library(kissDE)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)

list.files()

#Reading in data files from the Kissplice run
counts <- kissplice2counts("results_Huhu-S1_R1_Huhu-S1_R2_Huhu-S2_R1_Huhu-S2_R2_Huhu-S3_R1_Huhu-S3_R2_Huhu-L1_R1_Huhu-L1_R2_Huhu-L2_R1_Huhu-L2_R2_Huhu-L3_R1_Huhu-L3_R2_k41_coherents_type_0a.fa", counts = 0, pairedEnd = TRUE)
conditions <- c(rep("Huhu-S", 3), rep("Huhu-L", 3))
tiff('AS_counts_PCA.tiff', units="in", width=7.5, height=7.5, res=300, compression = 'lzw')
qualityControl(counts, conditions)
dev.off()
names(counts) 
str(counts)
head(counts$countsEvents)

#Reading in data files from the Kissplice run
results <- diffExpressedVariants(counts, conditions, nbCore = 3, pvalue = 0.05, filterLowCountsVariants = 5)
writeOutputKissDE(results, output = "kissDE_output.tsv")
exploreResults("kissDE_output.tsv.rds")

# Load KissDE results
kiss_de_results <- fread('kissDE_output.tsv')

# Inspect the first few rows of the dataset
head(kiss_de_results)

# Replace zero p-values with a small value
kiss_de_results$`15.Adjusted_pvalue`[kiss_de_results$`15.Adjusted_pvalue` == 0] <- .Machine$double.xmin

# Calculate -log10(Adjusted p-value) and remove points above 12
kiss_de_results <- kiss_de_results %>%
  mutate(log10_adj_pvalue = -log10(`15.Adjusted_pvalue`)) %>%
  filter(log10_adj_pvalue <= 12)

# Filter significant events
significant_events <- kiss_de_results %>%
  filter(`15.Adjusted_pvalue` < 0.05)

# Summarize significant events
kiss_de_summary <- significant_events %>%
  summarize(
    total_differential_events = n(),
    mean_DeltaPSI = mean(`16.Deltaf/DeltaPSI`),
    mean_adjusted_pvalue = mean(`15.Adjusted_pvalue`)
  )
print(kiss_de_summary)

# Identify top 10 significant points based on absolute Delta PSI
top10_significant <- significant_events %>%
  arrange(desc(abs(`16.Deltaf/DeltaPSI`))) %>%
  head(10)

# Add a column to classify points as upregulated or downregulated
kiss_de_results <- kiss_de_results %>%
  mutate(regulation = ifelse(`16.Deltaf/DeltaPSI` > 0, "Upregulated", "Downregulated"))

# Ensure there are no NA values in the regulation column
kiss_de_results <- kiss_de_results %>%
  mutate(regulation = ifelse(is.na(regulation), "Unclassified", regulation))

# Volcano plot
p_volcano <- ggplot(kiss_de_results, aes(x = `16.Deltaf/DeltaPSI`, y = log10_adj_pvalue, color = regulation)) +
  geom_point(alpha = 0.5, size = 2) +  # Increase the size of the dots
  geom_point(data = top10_significant, aes(x = `16.Deltaf/DeltaPSI`, y = log10_adj_pvalue), color = "black", size = 2) +  # Highlight top 10 points
  geom_text_repel(data = top10_significant, aes(label = `#1.ID`), color = "black", size = 3, box.padding = 0.5, max.overlaps = Inf) +  # Make labels bold and increase max overlaps
  labs(title = "Differential splicing events (Small vs Large)", x = "Delta PSI", y = "-log10(Adjusted p-value)", color = "Regulation") +
  scale_color_manual(values = c("Upregulated" = "goldenrod2", "Downregulated" = "gold2")) +  # Custom colors
  scale_x_continuous(name = expression(Delta~"PSI"), limits = c(-1, 1), 
                     breaks = seq(floor(min(kiss_de_results$`16.Deltaf/DeltaPSI`, na.rm = TRUE)), 
                                  ceiling(max(kiss_de_results$`16.Deltaf/DeltaPSI`, na.rm = TRUE)), by = 0.1)) +
  scale_y_continuous(name = expression(-log[10]~"Adjusted p-value"), limits = c(1.5, 12), 
                     breaks = seq(0, 13, by = 0.5)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    panel.grid.major = element_line(color = "grey90", size = 0.2),  # Lighter grid lines
    panel.grid.minor = element_line(color = "grey95", size = 0.1),  # Even lighter minor grid lines
    axis.line = element_line(color = "black", size = 0.5),  # Darker axis lines
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"))
print(p_volcano)

# Save the volcano plot
ggsave("volcano_plot_kissDE.tiff", plot = p_volcano, width = 7.5, height = 5.1, units = "in", dpi = 1200, compression = 'lzw')
