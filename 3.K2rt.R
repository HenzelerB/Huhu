library(tidyverse)
library(pheatmap)

# Reading in data files from the KissDE run
Huhu_k2rt <- read.table(file = 'Huhu_k2rt.tsv', sep = '\t', header = TRUE, check.names = FALSE)

# Removing unwanted annotations
Huhu_k2rt_filtered <- Huhu_k2rt %>%
  filter(
    Possible_sequencing_error %in% c(FALSE, "FALSE", "False", 0),
    Amino_acid_1 != 'N/A',  # Assuming columns for amino acids exist
    Amino_acid_2 != 'N/A',
    Is_condition_specific %in% c(TRUE, "TRUE", "True", 1),
    `KissDE_p-value` < 0.05,
    `KissDE_p-value` != 0,
    !KissDE_DeltaF %in% c("NA", "NaN", "N/A"),
    SNP_in_mutliple_assembled_genes %in% c(TRUE, "TRUE", "True", 1)
  )
head(Huhu_k2rt_filtered)


# Create a summary of amino acid changes
amino_acid_changes <- Huhu_k2rt_filtered %>%
  select(Amino_acid_1, Amino_acid_2) %>%
  group_by(Amino_acid_1, Amino_acid_2) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# Transform data for heat map
amino_acid_matrix <- amino_acid_changes %>%
  spread(key = Amino_acid_2, value = Count, fill = 0)

# Convert data to matrix
amino_acid_matrix_data <- as.matrix(amino_acid_matrix[,-1])
rownames(amino_acid_matrix_data) <- amino_acid_matrix$Amino_acid_1
colnames(amino_acid_matrix_data) <- colnames(amino_acid_matrix)[-1]

# Customize heat map colors and annotations
heatmap_colors <- colorRampPalette(c("ivory1", "#90EE90"))(256)

# Plotting the heat map with pheatmap
AA_change <- pheatmap(
  amino_acid_matrix_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = heatmap_colors,
  fontsize = 7,
  fontsize_row = 7,
  fontsize_col = 7,
  main = "Amino acid changes",
  display_numbers = FALSE,
  number_color = "black",
  border_color = "grey90",
  legend = TRUE,  # Ensure legend is displayed
  legend_breaks = c(0, 10, 20, 30, 40, 50),
  legend_labels = c("0", "10", "20", "30", "40", "50")
)

# Display the plot in the console
print(AA_change)

# Save the plot to a file
tiff('Aminoacid_change.tiff', units = "in", width = 4.5, height = 4.5, res = 1200, compression = 'lzw')
print(AA_change)
dev.off()
