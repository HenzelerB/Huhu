#Load libraries
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(readxl)

#Read the DESeq2 data
deseq2_data <- read.csv("DESeq2.csv")

#Extract species information from the column containing ID and group them
deseq2_data$Species <- str_extract(deseq2_data$ID, "\\[(.*?)\\]") %>% str_remove_all("[\\[\\]]")
deseq2_data <- deseq2_data %>%
  mutate(Regulation = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))

#Count the number of up- and down-regulated genes per species
regulation_counts <- deseq2_data %>%
  group_by(Species, Regulation) %>%
  summarise(count = n()) %>%
  spread(Regulation, count, fill = 0) %>%
  mutate(Ratio = `Up-regulated` / `Down-regulated`) %>%
  arrange(desc(`Up-regulated`), desc(`Down-regulated`)) %>%
  head(25) %>%
  column_to_rownames(var = "Species") %>%
  as.data.frame()
write.csv(regulation_counts, "species_DE.csv")

#Define color mapping for the `Down-regulated` and `Up-regulated` columns
col.s <- colorRamp2(c(0, max(regulation_counts, na.rm = TRUE) / 4, 
                      max(regulation_counts, na.rm = TRUE) / 3, 
                      max(regulation_counts, na.rm = TRUE) / 2, 
                      max(regulation_counts, na.rm = TRUE)),
                    c("ivory", "lightgoldenrod2", "springgreen3", "palegreen3", "salmon"))

#Define a separate color mapping for the `Ratio` column
ratio_col.s <- colorRamp2(c(0, max(regulation_counts$Ratio, na.rm = TRUE)), 
                          c("ivory", "goldenrod1"))

#Create the heatmap with `Ratio` as an additional column using HeatmapAnnotation
svg("DE_Species_Hmap.svg", width = 4, height = 8, bg = "transparent", family = "arial")

#Main heatmap for Down-regulated and Up-regulated columns with values displayed in cells
ht_main <- Heatmap(as.matrix(regulation_counts[, c("Down-regulated", "Up-regulated")]),
                   name = "Counts",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   col = col.s,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_names_gp = gpar(fontsize = 10),
                   row_names_gp = gpar(fontsize = 10),
                   row_names_side = "left",
                   column_names_rot = 90,
                   column_title = "Species distribution - Differential expression",
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                   heatmap_legend_param = list(title = "Counts"),
                   rect_gp = gpar(col = "grey", lwd = 0.5),
                   border_gp = gpar(col = "black", lty = 0.5),
                   # Adding values for Down-regulated and Up-regulated in cells
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (j == 1) {
                       grid.text(sprintf("%d", regulation_counts$`Down-regulated`[i]), x, y, gp = gpar(fontsize = 8, col = "black"))
                     } else if (j == 2) {
                       grid.text(sprintf("%d", regulation_counts$`Up-regulated`[i]), x, y, gp = gpar(fontsize = 8, col = "black"))
                     }
                   })

#Adding the Ratio column as a separate heatmap with values displayed in cells
ht_ratio <- Heatmap(as.matrix(regulation_counts[, "Ratio", drop = FALSE]),
                    name = "Ratio",
                    cluster_rows = FALSE,
                    col = ratio_col.s,
                    show_row_names = FALSE,
                    show_column_names = TRUE,
                    column_names_gp = gpar(fontsize = 10),
                    heatmap_legend_param = list(title = "Ratio"),
                    rect_gp = gpar(col = "grey", lwd = 0.5),
                    border_gp = gpar(col = "black", lty = 0.5),
                    # Adding the ratio values as text in cells
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", regulation_counts$Ratio[i]), x, y, gp = gpar(fontsize = 8, col = "black"))
                    })

#Draw the heatmaps together
draw(ht_main + ht_ratio)
dev.off()

#Read the datasets
df1 <- read_excel("Group_HS1-VS-HL1_DE_significant_anno.xlsx")
df2 <- read_excel("Group_HS2-VS-HL2_DE_significant_anno.xlsx")
df3 <- read_excel("Group_HS3-VS-HL3_DE_significant_anno.xlsx")

#Function to preprocess individual datasets
preprocess_dataset <- function(data) 
{
  data <- data %>%
    mutate(Regulation = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
    select(logFC, FDR, NR_Description, Regulation)
  return(data)
}

#Preprocess each dataset
df1 <- preprocess_dataset(df1)
df2 <- preprocess_dataset(df2)
df3 <- preprocess_dataset(df3)

# Merge datasets and calculate averages
m_df <- bind_rows(df1 %>% mutate(Dataset = "df1"), df2 %>% mutate(Dataset = "df2"), df1 %>% mutate(Dataset = "df3"))

#Calculate averages for logFC and FDR
avg_df <- m_df %>%
  group_by(NR_Description) %>%
  summarize(logFC = mean(logFC, na.rm = TRUE), FDR = mean(FDR, na.rm = TRUE), .groups = 'drop') %>%
  mutate(logFC = round(logFC, 4), Regulation = ifelse(logFC > 0, "Upregulated", "Downregulated"), CommonName = stringr::str_extract(NR_Description, "(?<=\\[).*?(?=\\])"), NR_Description = stringr::str_remove(NR_Description, "\\[.*?\\]"))

#Split the data into upregulated and downregulated groups
up <- avg_df %>%
  filter(Regulation == "Upregulated") %>%
  group_by(CommonName) %>%
  summarize(
    Count = n(),
    Avg_logFC = mean(logFC, na.rm = TRUE),
    .groups = 'drop')

down <- avg_df %>%
  filter(Regulation == "Downregulated") %>%
  group_by(CommonName) %>%
  summarize(
    Count = n(),
    Avg_logFC = mean(logFC, na.rm = TRUE),
    .groups = 'drop')

#Merge upregulated and downregulated groups by CommonName
Species_df <- full_join(up %>% rename(Upregulated_Count = Count, Upregulated_Avg_logFC = Avg_logFC),
                        down %>% rename(Downregulated_Count = Count, Downregulated_Avg_logFC = Avg_logFC),
                        by = "CommonName")%>%
  select(CommonName, Upregulated_Count, Downregulated_Count)
write.csv(Species_df,"Species_DE.csv")

#Load necessary libraries
library(tidyverse)
library(taxize)

#Function to fetch taxonomy for a species
fetch_taxonomy <- function(species_name) 
{tryCatch({result <- classification(species_name, db = "itis")
if (!is.null(result[[1]])) {
  taxonomy <- result[[1]]
  taxonomy$species <- species_name
  return(taxonomy)
} 
else 
{
  return(data.frame(rank = NA, name = NA, species = species_name))
}}, error = function(e) {
  return(data.frame(rank = NA, name = NA, species = species_name))})}

#Read the input file
input_file <- "DE.txt"
species_list <- readLines(input_file)

#Standardize taxonomy results for consistent columns
standardize_taxonomy <- function(taxonomy_df, species_name) {
  if (!is.data.frame(taxonomy_df)) {taxonomy_df <- as.data.frame(taxonomy_df)}
  taxonomy_df$species <- species_name
  all_columns <- c("rank", "name", "species")
  missing_columns <- setdiff(all_columns, colnames(taxonomy_df))
  for (col in missing_columns) 
  {
    taxonomy_df[[col]] <- NA
  }
  taxonomy_df <- taxonomy_df[, all_columns]
  return(taxonomy_df)}

#Combine taxonomy results into a single data frame
taxonomy_df <- do.call(rbind, lapply(names(taxonomy_results), function(species) 
{result <- taxonomy_results[[species]]
if (!is.null(result)) {
  return(standardize_taxonomy(result, species))
} 
else 
{
  return(data.frame(rank = NA, name = NA, species = species))}}))

#Save the taxonomy to a CSV file
output_file <- "taxonomy_results.csv"
write.csv(taxonomy_df, output_file, row.names = FALSE)

#Load necessary libraries
library(tidyverse)
library(circlize)
library(pheatmap)
library(ComplexHeatmap)

#Load the datasets
species_de <- read.csv("Species_DE.csv")
class_de <- read.csv("class_DE.csv")

#Merge the datasets
merged_df <- species_de %>%
  inner_join(class_de, by = c("CommonName" = "Species")) %>%
  select(-contains("Unnamed")) %>%
  group_by(Name) %>%
  summarise(Sp.Count = n(), 'Up-reg' = as.integer(sum(Upregulated_Count, na.rm = TRUE)), 'Down-reg' = as.integer(sum(Downregulated_Count, na.rm = TRUE))) %>%
  arrange(desc(Sp.Count)) %>%
  select('Name', 'Down-reg', 'Up-reg', 'Sp.Count')%>%
  column_to_rownames("Name") %>%
  head(25)%>%
  as.matrix()

#Color coding for the heatmap
col.s <- colorRamp2(c(0, max(merged_df, na.rm = TRUE) / 5, max(merged_df, na.rm = TRUE) / 4, max(merged_df, na.rm = TRUE) / 3), c("ivory", "goldenrod1", "palegreen3", "salmon"))
svg("DE_Class.svg", width = 4, height = 8, bg = "transparent", family = "arial")
Heatmap(merged_df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col.s,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left",
        column_names_rot = 90,  
        column_title = "Differentially expressed classes",  
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(title = "Counts"),
        rect_gp = gpar(col = "grey", lwd = 0.5),
        border_gp = gpar(col = "black", lty = 0.5),
        column_gap = unit(2, "mm"),
        column_title_rot = 0,
        # Add cell values as text within each heatmap block
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", merged_df[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        })
dev.off()
