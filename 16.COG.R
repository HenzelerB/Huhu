#Load necessary libraries
library(readxl)
library(dplyr)
library(stringr)

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
write.csv(m.df,"Functional_Annotations.csv")

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
  select(-Gene_ID, -Original_Uniprot_ID)%>%
  rename(Gene_ID = Entry.Name)%>%
  rename(Gene_Name = Gene.Names.Primary)%>%
  rename(S.No = X)%>%
  select(S.No, Gene_ID, Gene_Name, Description, length, Species, logFC, logCPM, pval, FDR, Regulation, GO_BP, GO_CC, GO_MF, KO_ID, KO_Definition, COG_ID, COG_Description)
write.csv(merged_data,"Annotations.csv")

#Separate COG Analysis dataset
cog <- merged_data %>%
  select(Gene_ID, Gene_Name, COG_ID, COG_Description, logFC, logCPM, pval, FDR, Regulation) %>%
  filter(!is.na(COG_ID))
write.csv(cog, "COG.csv", row.names = FALSE)
