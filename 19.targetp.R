# Load necessary libraries
library(tidyverse)
library(VennDiagram)
library(readxl)
library(svglite)

# Read and preprocess the annotation data
Ann <- read_excel("Annotation.xlsx") %>%
  select(geneID, Uniprot_Acc, Uniprot_Description)

# Process the Small dataset
column_names <- c("ID", "Prediction", "SP", "noTP", "mTP", "CS_Position")
df.S <- read.table('Small_P_output_summary.targetp2', 
                   sep = "\t",
                   header = FALSE,
                   col.names = column_names,
                   fill = TRUE,
                   quote = "",
                   stringsAsFactors = FALSE) %>%
  mutate(ID = str_extract(ID, "TRINITY_[^\\s]+") %>% str_remove("\\.p[1-5]$")) %>%
  mutate(geneID= gsub("^TRINITY_", "", ID)) %>%
  left_join(Ann, by = c("geneID" = "geneID")) %>%
  filter(!is.na(Uniprot_Acc) & !Uniprot_Acc %in% c("-", "--", "0"))%>%
  select(-ID)
df.L <- read.table('Large_P_output_summary.targetp2', 
                   sep = "\t",
                   header = FALSE,
                   col.names = column_names,
                   fill = TRUE,
                   quote = "",
                   stringsAsFactors = FALSE) %>%
  mutate(ID = str_extract(ID, "TRINITY_[^\\s]+") %>% str_remove("\\.p[1-5]$")) %>%
  mutate(geneID= gsub("^TRINITY_", "", ID)) %>%
  left_join(Ann, by = c("geneID" = "geneID")) %>%
  filter(!is.na(Uniprot_Acc) & !Uniprot_Acc %in% c("-", "--", "0"))%>%
  select(-ID)

# Extract gene IDs from both datasets
large_geneIDs <- unique(df.L$geneID)
small_geneIDs <- unique(df.S$geneID)

# Create a Venn diagram
pdf("targetP.pdf", width = 10, height = 10, onefile = TRUE)
venn.plot <- draw.pairwise.venn(area1 = length(large_geneIDs),
                                area2 = length(small_geneIDs),
                                cross.area = length(intersect(large_geneIDs, small_geneIDs)),
                                category = c("Large Dataset", "Small Dataset"),
                                fill = c("goldenrod1", "lightgoldenrod"),
                                alpha = 0.5,
                                cat.pos = c(0, 0),
                                cat.dist = 0.03)
dev.off()

# Filter geneIDs with mTP > 0.8
L_genes_above_0_8 <- df.L %>%
  filter(mTP > 0.8) %>%
  select(Uniprot_Description, mTP, SP, noTP, CS_Position) %>%
  distinct()
S_genes_above_0_8 <- df.S %>%
  filter(mTP > 0.8) %>%
  select(Uniprot_Description, mTP, SP, noTP, CS_Position) %>%
  distinct()

# Function to calculate counts
calculate_counts <- function(data, cutoff = 0.8) {
  Secretory <- sum(data$SP > cutoff, na.rm = TRUE)
  `Cytosolic/non-organelle protein` <- sum(data$noTP > cutoff, na.rm = TRUE)
  Mitochondrial <- sum(data$mTP > cutoff, na.rm = TRUE)
  return(data.frame(
    `Secretory` = Secretory, 
    `Cytosolic/non-organelle protein` = `Cytosolic/non-organelle protein`, 
    `Mitochondrial` = Mitochondrial
  ))
}

# Add a column to indicate the dataset source
"Large" <- L_genes_above_0_8$Dataset%>%
  head(15)
"Small" <- S_genes_above_0_8$Dataset%>%
  head(15)

# Combine the datasets
combined_genes <- rbind(L_genes_above_0_8%>%head(15), S_genes_above_0_8%>%head(15))

# Select relevant columns (include Dataset column)
plot_data <- combined_genes[, c("Uniprot_Description", "mTP", "Dataset", "SP")]

# Order data by Uniprot_Description
plot_data <- plot_data[order(plot_data$Dataset, plot_data$Uniprot_Description), ]

# Create a rank column for y-axis
plot_data$Rank <- seq_along(plot_data$mTP)

# Create the plot with size based on SP/Secretory value
pdf("mTP.pdf", width = 10, height = 6, onefile = TRUE)
ggplot(plot_data, aes(x = mTP, y = Rank, color = Dataset, size = SP)) +
  geom_point(alpha = 0.7) +  # Add transparency to better visualize overlapping points
  scale_y_reverse(breaks = plot_data$Rank, labels = plot_data$Uniprot_Description) +
  facet_wrap(~Dataset, scales = "free_y", ncol = 1) + # Separate panels for each dataset
  scale_color_manual(
    values = c("Large" = "goldenrod1", "Small" = "lightgoldenrod") # Define custom colors
  ) +
  scale_size_continuous(
    name = "SP (Secretory)", 
    range = c(1, 5)  # Adjust the size range for better visualization
  ) +
  scale_x_continuous(
    breaks = seq(floor(min(plot_data$mTP)), ceiling(max(plot_data$mTP)), by = 0.025), # Adjust the step size
    minor_breaks = NULL # Optional: Remove minor breaks for a cleaner look
  ) +
  labs(
    title = "Uniprot Descriptions and mTP Scores (Grouped by Dataset)",
    x = "mTP Score",
    y = "Proteins (ranked by mTP score)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray", linetype = "dashed", size = 0.05),
    strip.text = element_text(size = 12, face = "bold") # Format facet titles
  )
dev.off()

