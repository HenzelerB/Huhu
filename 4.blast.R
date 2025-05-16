# Load necessary libraries
library(tidyverse)
library(data.table)
library(cluster)
library(writexl)
library(ggdendro)
library(extrafont)

# Read the data
Huhu_blastn <- fread('Nr_Annotation.xls', sep = '\t', header = FALSE, fill = TRUE)
colnames(Huhu_blastn) <- c("Gene_ID", "NR_protein_ID", "Description", "Identity", "Evalue", "Score")

# Extract and clean the Species_Name
Huhu_blastn <- Huhu_blastn %>%
  mutate(Species_Name = str_extract(Description, "\\[([^\\]]+)\\]") %>% str_remove_all("\\[|\\]")) %>%
  mutate(Species_Name = str_replace(Species_Name, "(\\w+)\\s+(\\w+)", "\\1.\\2") %>%
           str_replace("^(\\w)\\w*\\.", "\\1.") %>%
           str_to_title())

# Save the updated data to a CSV file
fwrite(Huhu_blastn, 'Updated_Nr_Annotation.csv')

# Calculate species counts and select the top 15 species
species_counts <- Huhu_blastn %>%
  group_by(Species_Name) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

top_15_species <- species_counts %>%
  slice_head(n = 15)

# Filter the original data to include only the top 15 species
Huhu_blastn_top_15 <- Huhu_blastn %>%
  filter(Species_Name %in% top_15_species$Species_Name)

# Ensure Identity is numeric
Huhu_blastn_top_15 <- Huhu_blastn_top_15 %>%
  mutate(Identity = as.numeric(str_remove(Identity, "%")))

# Aggregating data by species and calculating mean Identity for the top 15 species
aggregated_data <- Huhu_blastn_top_15 %>%
  group_by(Species_Name) %>%
  summarise(mean_identity = mean(Identity, na.rm = TRUE)) %>%
  filter(!is.na(mean_identity))

# Clustering the species based on the mean Identity
clustering <- hclust(dist(aggregated_data$mean_identity))

# Export to PDF file with Arial font
svg('Cluster.svg', width = 8, height = 8, family = "Arial")  # Specify font family as Arial
plot(clustering, labels = aggregated_data$Species_Name, main = "", xlab = "", sub = "", ylab = "", yaxt = "n")

# Add distance labels to the dendrogram branches
for (i in 1:length(clustering$height)) {
  x_center <- mean(clustering$merge[i, ])
  y_height <- clustering$height[i]
  text(x_center, y_height, round(y_height, 1), pos = 3, cex = 0.8, family = "Arial")
}

# Highlight clusters with k = 3
rect.hclust(clustering, k = 3, border = c("grey", "grey", "goldenrod1"))

# Close the device
dev.off()
