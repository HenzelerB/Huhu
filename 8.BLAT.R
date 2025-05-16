# Load necessary libraries
library(tidyverse)
library(data.table)
library(ggtext)  # For element_markdown

# Define column names for PSL format
psl_colnames <- c("matches", "misMatches", "repMatch", "nCount", "qNumInsert", 
                  "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", 
                  "qName", "qSize", "qStart", "qEnd", "tName", "tSize", 
                  "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")

# Read the PSL file using fread with the correct delimiter and column names
BLAT <- fread(file = 'Huhu.psl', sep = '\t', header = FALSE, col.names = psl_colnames)

# Filter out the 'qName' column
BLAT_filtered <- BLAT %>% select(qStart, qEnd)

# Define custom breaks for the axes
x_breaks <- seq(0, 100, by = 10)
y_breaks <- seq(0, 60, by = 10)

# Create a scatter plot with reversed axes (qEnd vs qStart) and customized axis units
qEnd_qStart_plot <- ggplot(BLAT_filtered, aes(x = qEnd, y = qStart)) +
  geom_point(alpha = 0.5, color = "goldenrod1", size = 0.5) +
  labs(title = "qEnd vs qStart",
       x = "qEnd",
       y = "qStart") +
  theme_minimal() +
  scale_x_continuous(breaks = x_breaks, limits = c(20, 90)) +
  scale_y_continuous(breaks = y_breaks, limits = c(0, 60))

# Save the plot as a TIFF image
tiff('qEnd_qStart_plot.tiff', units="in", width=8, height=3, res=1200, compression = 'lzw')
print(qEnd_qStart_plot)
dev.off()
