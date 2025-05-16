#Load required libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(svglite)

#Load the datasets
go <- read.csv("enriched_go_uniprot.csv")%>%
  filter(!(GO_Term %in% c("-", "--", "0", "")) & !is.na(GO_Term)) %>%
  filter(rowSums(is.na(.) | . == "") < ncol(.))

#Function to split and restructure data for a GO category
df <- go %>%
  group_by(GO_Term) %>%
  summarise(GO_Category = first(GO_Category), Protein_Names = first(Protein.names), Entry_Names = paste(unique(Entry.Name), collapse = ","), Gene_Counts = length(unique(Entry.Name)))%>%
  arrange(desc(Gene_Counts))

#Process sum of GO category
category_counts <- df %>%
  group_by(GO_Category) %>%
  summarise(Counts = n()) %>%
  arrange(desc(Counts)) %>%
  mutate(Percentage = round((Counts / sum(Counts)) * 100, 1), GO_Category = recode(GO_Category, "Biological Process" = "BP", "Cellular Component" = "CC", "Molecular Function" = "MF"))

#Create a pie chart
custom_colors <- c("BP" = "salmon", "CC" = "goldenrod1", "MF" = "springgreen1")
pie_chart <- ggplot(category_counts, aes(x = "", y = Counts, fill = GO_Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = custom_colors) +
  geom_text(aes(label = paste0(Counts, " ")), 
            position = position_stack(vjust = 0.5)) +
  labs(title = "GO distribution", x = NULL, y = NULL, fill = "Category") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("go_category.svg", plot = pie_chart, width = 5, height = 5)

#Define a function to process datasets
process_dataset <- function(file_path) 
{read_excel(file_path) %>%
    select(Gene_ID, length, logFC, logCPM, pval, FDR, Uniprot_Acc, NR_Description, Regulation, GO_BP, GO_CC, GO_MF, KO_ID, KO_Definition, COG_ID, COG_Description) %>%
    mutate(Species = str_extract(NR_Description, "\\[.*?\\]"), Species = str_remove_all(Species, "\\[|\\]"), NR_Description = str_remove(NR_Description, "\\[.*?\\]")) %>%
    mutate(Species = str_trim(Species), Description = str_trim(NR_Description)) %>%
    select(-NR_Description) %>%
    filter(if_all(everything(), ~ !is.na(.) & . != "-" & . != "--" & . != "0" & . != ""))}

#Process each dataset using the function
df1 <- process_dataset("Group_HS1-VS-HL1_DE_significant_anno.xlsx")
df2 <- process_dataset("Group_HS2-VS-HL2_DE_significant_anno.xlsx")
df3 <- process_dataset("Group_HS3-VS-HL3_DE_significant_anno.xlsx") %>% rename_with(~ paste0(., "_3"), -Gene_ID)

#Merge the three datasets based on Gene_ID
m.df <- df1 %>%
  full_join(df2, by = "Gene_ID", suffix = c("_1", "_2")) %>%
  full_join(df3, by = "Gene_ID") %>%
  rowwise() %>%
  mutate(
    logFC_avg = mean(c_across(matches("^logFC_")), na.rm = TRUE),
    logCPM_avg = mean(c_across(matches("^logCPM_")), na.rm = TRUE),
    pval_avg = mean(c_across(matches("^pval_")), na.rm = TRUE),
    FDR_avg = mean(c_across(matches("^FDR_")), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(Uniprot_Acc_1) & Uniprot_Acc_1 != "0" & Uniprot_Acc_1 != "-" & Uniprot_Acc_1 != "--") %>%
  select(Gene_ID, logFC_avg, logCPM_avg, pval_avg, FDR_avg, Uniprot_Acc_1)

#Rejoin NR_Description and Species columns from the initial df1
m.df <- m.df %>%
  left_join(df1 %>% select(Gene_ID, Description, Species, length, Regulation, GO_BP, GO_CC, GO_MF, KO_ID, KO_Definition, COG_ID, COG_Description), by = "Gene_ID")%>%
  rename(Uniprot_ID = Uniprot_Acc_1)%>%
  rename(logFC = logFC_avg)%>%
  rename(logCPM = logCPM_avg)%>%
  rename(pval = pval_avg)%>%
  rename(FDR = FDR_avg)

#Load the data
functional_annotations <- read.csv("Functional_Annotations.csv", stringsAsFactors = FALSE)
id_mapping <- read.delim("idmapping_2024_11_28.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Rename columns for clarity
colnames(id_mapping) <- c("From", "Entry", "Entry.Name", "Protein.Names", 
                          "Gene.Names", "Gene.Names.Synonym", "Gene.Names.Primary")
#Merge datasets, keeping both identifiers
merged_data <- functional_annotations %>%
  rename(Original_Uniprot_ID = Uniprot_ID)%>%
  left_join(id_mapping %>% select(Entry, Entry.Name, Gene.Names.Primary), by = c("Original_Uniprot_ID" = "Entry"))%>%
  select(-Gene_ID)%>%
  rename(Gene_ID = Entry.Name)%>%
  rename(Uniprot_ID = Original_Uniprot_ID)%>%
  rename(Gene_Name = Gene.Names.Primary)%>%
  rename(S.No = X)%>%
  select(S.No, Gene_ID, Uniprot_ID, Gene_Name, Description, length, Species, logFC, logCPM, pval, FDR, Regulation, GO_BP, GO_CC, GO_MF, KO_ID, KO_Definition, COG_ID, COG_Description)

#Separate GO Analysis datasets
go <- merged_data %>%
  select(Gene_ID, Uniprot_ID, Gene_Name, GO_BP, GO_CC, GO_MF, logFC, logCPM, pval, FDR, Regulation) %>%
  filter(!is.na(GO_BP) & !is.na(GO_CC) & !is.na(GO_MF))

#Function to process GO terms by regulation status with statistics
  process_go_terms <- function(go, go_column) {
    processed_data <- go %>%
      select(!!sym(go_column), Gene_ID, Uniprot_ID, Gene_Name, logFC, pval, logCPM, FDR) %>%
      filter(!is.na(!!sym(go_column))) %>%
      separate_rows(!!sym(go_column), sep = ";") %>%
      group_by(!!sym(go_column)) %>%
      summarise(
        Gene_ID = paste(unique(Gene_ID), collapse = ","),
        Gene_Name = paste(unique(Gene_Name), collapse = ","),
        Uniprot_ID = paste(unique(Uniprot_ID), collapse = ","),
        .groups = "drop"
      )
    colnames(processed_data)[1] <- "GO_Term"
    return(processed_data)
  }

#Process each GO category by regulation with statistics
go_bp <- process_go_terms(go, "GO_BP") %>%
  mutate(GO_Category = "Biological Process")
go_cc <- process_go_terms(go, "GO_CC") %>%
  mutate(GO_Category = "Cellular Component")
go_mf <- process_go_terms(go, "GO_MF") %>%
  mutate(GO_Category = "Molecular Function")

#Combine all GO terms into a single dataframe
go_terms <- bind_rows(go_bp, go_cc, go_mf)%>%
  mutate(Gene_Count = sapply(strsplit(Gene_ID, ","), length))%>%
  filter(!is.na(GO_Term) & GO_Term != "" & GO_Term != "-" & GO_Term != "--")%>%
  mutate(GO_Category = sub("(^GO:\\d+):.*", "\\1", GO_Term), Term = sub("^GO:\\d+:(.*)", "\\1", GO_Term))%>%
  select(-GO_Term)%>%
  select(Gene_ID, Gene_Name, Uniprot_ID,Gene_Count, GO_Category, Term)
write.csv(go_terms, "go_annotations.csv")

#Read the datasets
df1 <- read_excel("HS1-VS-HL1.xlsx")
df2 <- read_excel("HS2-VS-HL2.xlsx")
df3 <- read_excel("HS3-VS-HL3.xlsx")

#Define common columns for merging
common_columns <- c("category", "term", "ontology")

#merge the datasets and convert numeric columns to numeric type
m.df <- df1%>%
  inner_join(df2, by = common_columns, suffix = c("_1", "_2"))%>%
  inner_join(df3%>% rename_with(~ paste0(., "_3"), -common_columns), by = common_columns)%>%
  mutate(across(ends_with("_1"), as.numeric, .names = "{.col}"))%>%
  mutate(across(ends_with("_2"), as.numeric, .names = "{.col}"))%>%
  mutate(across(ends_with("_3"), as.numeric, .names = "{.col}"))

#Compute averages for each group of columns with matching suffixes
final_df <- m.df %>%
  mutate(numDEInCat_average = rowMeans(select(., numDEInCat_1, numDEInCat_2, numDEInCat_3), na.rm = TRUE),
         numInCat_average = rowMeans(select(., numInCat_1, numInCat_2, numInCat_3), na.rm = TRUE),
         over_represented_pvalue_average = rowMeans(select(., over_represented_pvalue_1, over_represented_pvalue_2, over_represented_pvalue_3), na.rm = TRUE),
         over_represented_FDR_average = rowMeans(select(., over_represented_FDR_1, over_represented_FDR_2, over_represented_FDR_3), na.rm = TRUE),
         GeneNumber_Up_average = rowMeans(select(., `GeneNumber(Up)_1`, `GeneNumber(Up)_2`, `GeneNumber(Up)_3`), na.rm = TRUE),
         GeneNumber_Down_average = rowMeans(select(., `GeneNumber(Down)_1`, `GeneNumber(Down)_2`, `GeneNumber(Down)_3`), na.rm = TRUE)) %>%
  select(all_of(common_columns), ends_with("_average"))
write.csv(final_df, "GO_datasets_merged.csv")

#Read the datasets
go_annotations <- read.csv("go_annotations.csv")
go_datasets_merged <- read.csv("GO_datasets_merged.csv")

#Merge the datasets by the common columns
merged_data <- merge(go_annotations, go_datasets_merged, by.x = "GO_Category", by.y = "category", all = FALSE)%>%
  select(Gene_ID, Gene_Name, Uniprot_ID, Gene_Count, GO_Category, 
         term, ontology, numDEInCat_average, numInCat_average, 
         over_represented_pvalue_average, over_represented_FDR_average, 
         GeneNumber_Up_average, GeneNumber_Down_average) %>%
  rename(numDEInCat = numDEInCat_average,
         numInCat = numInCat_average,
         over_represented_pvalue = over_represented_pvalue_average,
         over_represented_FDR = over_represented_FDR_average,
         GeneNumber_Up = GeneNumber_Up_average,
         GeneNumber_Down = GeneNumber_Down_average) %>%
  mutate(numDEInCat = round(numDEInCat),
         numInCat = round(numInCat),
         GeneNumber_Up = round(GeneNumber_Up),
         GeneNumber_Down = round(GeneNumber_Down),
         Ratio = numDEInCat / numInCat) %>%
  arrange(desc(Gene_Count))

#Save the merged data to a new CSV file
write.csv(merged_data, "merged_GO_datasets.csv", row.names = FALSE)

#Load the dataset
df <- read.csv("merged_GO_datasets.csv", stringsAsFactors = FALSE) 
df.go <- read.csv("go_terms_clusters.csv", stringsAsFactors = FALSE)
df.annotation <- df %>%
  filter(over_represented_FDR < 0.05, over_represented_pvalue < 0.05) %>%
  select(GO_Category, term, numDEInCat, over_represented_pvalue, ontology, GeneNumber_Up, GeneNumber_Down, Ratio) %>%
  rename(Code = GO_Category,
         Term = term,
         numDE = numDEInCat,
         pvalue = over_represented_pvalue,
         Type = ontology,
         Up = GeneNumber_Up,
         Down = GeneNumber_Down,
         Ratio = Ratio) %>%
  mutate(Type = recode(Type,
                       "BP" = "Biological process",
                       "CC" = "Cellular component",
                       "MF" = "Molecular function")) %>%
  arrange(desc(numDE)) %>%
  group_by(Type) %>%
  slice_max(numDE, n = 25) %>%
  ungroup() %>%
  arrange(desc(numDE))
write.csv(df.annotation,"GOPlot.csv")
go.tmp <- df.annotation %>%
  left_join(df.go %>% select(Term, Cluster, `Cluster.Name`), by = "Term")
write.csv(go.tmp, "go_terms_clusters.csv")
go <- read.csv("clusters.csv") %>%
  mutate(Term_With_Expression = paste0(Code, " (", Down, "↓)"," (", Up, "↑)"))

#Convert the long sheet to a matrix for calculating the distance.
go_matrix <- dcast(go,Term~Cluster, value.var = "numDE")
rownames(go_matrix) <- go_matrix[["Term"]]
go_matrix = go_matrix [,-1]
go_matrix[is.na(go_matrix)] = 0
Segement = data.frame()
Label = data.frame()
Num = 0
for( Level in unique(go$Type)){
  GROUP = go$Term[go$Type == Level]
  hc <- hclust(dist(go_matrix[row.names(go_matrix) %in% GROUP,]))
  dendr    <- dendro_data(hc, type="rectangle") # convert for ggplot
  S_tmp <- segment(dendr)
  S_tmp$Level <- Level
  L_tmp <- label(dendr)
  Segement <- rbind(Segement, S_tmp)
  Label <- rbind(Label, L_tmp)
}
colnames(Segement)[5] = "Type"

#Synchrony the Levels
go$Cluster = factor(go$Cluster, levels = c("hclust", as.character(unique(go$Cluster) )))
Segement$Cluster = "hclust"

#Normalised GOPlot
go.plot <- ggplot() +
  geom_point(data= go, aes(x = Cluster, y = Term, size = numDE, color = Ratio)) +
  facet_grid(Type ~ Cluster, scales = 'free', switch = "y") +
  scale_color_gradient2(high = "salmon", mid = "goldenrod1", low = "springgreen1", name = "Ratio") + 
  theme_bw() +   
  theme(
    strip.text = element_text(face = "bold"),   
    strip.background = element_rect(colour = "black", fill = FALSE),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text.y.right = element_text(hjust = 1)
  ) +
  labs(title = "GO Enrichment", x = "Clusters", y = "Categories", size = "Counts") +
  geom_segment(data = Segement, aes(x = -y, y = x, xend = -yend, yend = xend)) +
  scale_y_discrete(
    position = "right",
    labels = function(x) {
      # Map unique labels for each facet
      go %>%
        filter(Term %in% x) %>%
        arrange(match(Term, x)) %>%
        pull(Term_With_Expression)
    }
  )
ggsave(filename = "goplot.svg", device = "svg", plot = go.plot, scale = 2, height = 15, width = 10, units = "cm", dpi = 300)

