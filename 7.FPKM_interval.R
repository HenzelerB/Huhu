# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Read the file as a tab-separated values (TSV) file
FPKM <- read.table('interval.xls', sep = '\t', header = TRUE, fill = TRUE, quote = "")
head(FPKM)

# Melt the data for ggplot2
FPKM_melted <- gather(FPKM, key = "Interval", value = "Count", -FPKM.Interval)

# Convert the Count values from percentage format to numeric values
FPKM_melted <- FPKM_melted %>%
  mutate(Count = as.numeric(gsub("\\(.*\\)|%", "", Count)))
head(FPKM_melted)

# Remove the 'X' from the Interval values, replace '.' with spaces, '..' with '>', and 'to' with ' - '
FPKM_melted$Interval <- FPKM_melted$Interval %>%
  gsub("^X", "", .) %>%
  gsub("\\.\\.", " > ", .) %>%
  gsub("\\.", " ", .) %>%
  gsub(" to ", " - ", .)
head(FPKM_melted)

# Calculate total counts for each interval
interval_totals <- FPKM_melted %>%
  group_by(Interval) %>%
  summarise(Total = sum(Count)) %>%
  arrange(Total)

# Reorder intervals based on total counts
FPKM_melted$Interval <- factor(FPKM_melted$Interval, levels = interval_totals$Interval)

# Define pastel colors for fill
interval_colors <- c("#FFFF99", "#ADD8E6", "#90EE90", "#FFB6C1", "#FFC107")  # Corresponding to intervals

# Create the stacked bar plot with borders and transparency
FPKM <- ggplot(FPKM_melted, aes(fill = Interval, y = Count, x = FPKM.Interval)) +
  geom_bar(position = "stack", stat = "identity", width = 0.6, alpha = 1, colour = "grey60", size = 0.1) +  # Adjust border thickness here
  theme_minimal() +
  scale_fill_manual(values = interval_colors) +
  labs(title = "FPKM distribution",
       x = NULL,
       y = "Counts (%)") +
  theme(
    text = element_text(family = "Arial", size = 10),  # Set all text to Arial and size 10
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),
    axis.text.y = element_text(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),  # Title slightly larger
    legend.title = element_text(face = "bold"),  # Make legend title bold
    panel.grid.major = element_line(color = "grey70", size = 0.1),
    panel.grid.minor = element_line(color = "grey80", size = 0.1),
    panel.border = element_rect(color = "grey50", fill = NA, size = 0.2)
  )

print(FPKM)

# Save the plot to a file
ggsave("FPKM_Interval_Distribution_Ordered.tiff", plot = FPKM, width = 5, height = 3, dpi = 1200)
