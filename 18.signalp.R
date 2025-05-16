# Load necessary libraries
library(dplyr)
library(stringr)
library(readxl)
library(VennDiagram)
library(grid)
library(svglite)

# Read and preprocess the annotation data
Ann <- read_excel("Annotation.xlsx") %>%
  select(geneID, Uniprot_Acc, Uniprot_Description)

# Process the Small dataset
df.S <- read.csv("S.csv", header = TRUE) %>%
  mutate(ID = str_extract(X..ID, "TRINITY_[^\\s]+") %>% str_remove("\\.p[1-5]$")) %>%
  mutate(ID = gsub("^TRINITY_", "", ID)) %>%
  select(-X..ID) %>%
  left_join(Ann, by = c("ID" = "geneID")) %>%
  filter(!is.na(Uniprot_Acc) & !Uniprot_Acc %in% c("-", "--", "0"))
write.csv(df.S, "Small.csv")

# Process the Large dataset
df.L <- read.csv("L.csv", header = TRUE) %>%
  mutate(ID = str_extract(X..ID, "TRINITY_[^\\s]+") %>% str_remove("\\.p[1-5]$")) %>%
  mutate(ID = gsub("^TRINITY_", "", ID)) %>%
  select(-X..ID) %>%
  left_join(Ann, by = c("ID" = "geneID")) %>%
  filter(!is.na(Uniprot_Acc) & !Uniprot_Acc %in% c("-", "--", "0"))
write.csv(df.L, "Large.csv")

# Summarize a grouped dataset
summarize_grouped_data <- function(data, group_column) {
  data %>% group_by(!!sym(group_column)) %>% summarise(row_count = n())
}

# Summarize the datasets
summary.L <- summarize_grouped_data(df.L, "Prediction")
summary.S <- summarize_grouped_data(df.S, "Prediction")

# Identify common and unique Uniprot_Acc values
common <- intersect(df.L$Uniprot_Acc, df.S$Uniprot_Acc)
unique.L <- setdiff(df.L$Uniprot_Acc, df.S$Uniprot_Acc)
unique.S <- setdiff(df.S$Uniprot_Acc, df.L$Uniprot_Acc)

# Filter datasets to keep only unique rows
large_unique_df <- df.L[df.L$Uniprot_Acc %in% unique.L, ]
small_unique_df <- df.S[df.S$Uniprot_Acc %in% unique.S, ]

# Summarize the unique datasets
u.summary.L <- summarize_grouped_data(large_unique_df, "Prediction")
u.summary.S <- summarize_grouped_data(small_unique_df, "Prediction")

# Reload the datasets
large <- read.csv("Large.csv")
small <- read.csv("Small.csv")

# Create a function to generate a Venn diagram for a specific Prediction value
plot_venn_for_prediction <- function(prediction_type, large, small) {
  # Filter the datasets for the specific Prediction type
  large_filtered <- large %>% filter(Prediction == prediction_type) %>% pull(ID) %>% unique()
  small_filtered <- small %>% filter(Prediction == prediction_type) %>% pull(ID) %>% unique()
  
  # Save the Venn diagram as a PDF file
  file_name <- paste0(" ", prediction_type, ".pdf")
  pdf(file_name, width = 5, height = 5)
  
  # Create a new page for the grid
  grid.newpage()
  
  # Add a title
  grid.text(
    label = paste(" ", prediction_type),
    x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  # Create the Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(large_filtered),
    area2 = length(small_filtered),
    cross.area = length(intersect(large_filtered, small_filtered)),
    category = c("Large", "Small"),
    fill = c("darkgoldenrod1", "lightgoldenrod1"),
    cat.cex = 1.2,
    cex = 1.4
  )
  
  # Finalize and save the PDF
  dev.off()
  cat(paste(" ", prediction_type, "saved to", file_name, "\n"))
}

# Generate and save individual Venn diagrams for specific Prediction values
plot_venn_for_prediction("LIPO", large, small)
plot_venn_for_prediction("OTHER", large, small)
plot_venn_for_prediction("PILIN", large, small)
plot_venn_for_prediction("SP", large, small)
plot_venn_for_prediction("TAT", large, small)
plot_venn_for_prediction("TATLIPO", large, small)
