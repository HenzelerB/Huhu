#Loading libraries for the run
library(tidyr)
library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(ggsci)
library(gplots)
library(EnhancedVolcano)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(genefilter)
library(factoextra)
library(ggplot2)
library(ggrepel)  # For better text labels

# Read the datafile and adjust it
fcData = read.table('all.fpkm_anno.txt', sep='\t', header=TRUE)
fcData %>% head()
dim(fcData)
names(fcData)

# rename the data columns
names(fcData)[3:8] = c("Huhu-L1", "Huhu-L2", "Huhu-L3", "Huhu-S1", "Huhu-S2", "Huhu-S3")
fcData %>% head()

# Extract counts data
counts = fcData[, 3:8]
rownames(counts) = fcData$Gene_ID
counts %>% head()
counts_L = fcData[, 3:5]
rownames(counts_L) = fcData$Gene_ID
counts_L %>% head()
counts_S = fcData[, 6:8]
rownames(counts_S) = fcData$Gene_ID
counts_S %>% head()

#Read counts per sample
colSums(counts)
colSums(counts_L)
colSums(counts_S)

#Color coding the groups
Colour=c("goldenrod1", "goldenrod1", "goldenrod1", "goldenrod1", "goldenrod1", "goldenrod1")
Colour
Colour_small=c("goldenrod1", "goldenrod1", "goldenrod1")
Colour_small
Colour_large=c("lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1")
Colour_large

#Fixing the font size of
cex.main=2 #change font size of title
cex.sub=2 #change font size of subtitle
cex.lab=2 #change font size of axis labels
cex.axis=2 #change font size of axis text

#Visualizing the counts data
boxplot(as.matrix(counts) ~ col(counts), 
        key=F, trace="none", col=Colour, margin=c(5, 5), 
        xlab="Samples", 
        ylab="Counts", main="Counts for small and large", cex.main=0.8, cex.sub=0.8, cex.lab=0.8, cex.axis=0.8)

# Check the dimensions before removal
dim(counts)
dim(counts_L)
dim(counts_S)

#Normalizing the zero counts from the set
colSums(counts==0)
colSums(counts_L==0)
colSums(counts_S==0)

# Identify rows with zero counts
Counts <- rowSums(counts == 0) > 0
Counts_L <- rowSums(counts_L == 0) > 0
Counts_S <- rowSums(counts_S == 0) > 0

# Remove rows with zero counts
counts <- counts[!Counts, ]
counts_L <- counts_L[!Counts_L, ]
counts_S <- counts_S[!Counts_S, ]

# Check the dimensions after removal
dim(counts)
dim(counts_L)
dim(counts_S)

#Log transformation (add 0.5 to avoid log(0) issues)
logCounts = log2(as.matrix(counts)+ 1.5)
logCounts_L = log2(as.matrix(counts_L)+ 1.5)
logCounts_S = log2(as.matrix(counts_S)+ 1.5)

#Visualizing the log-transformed counts data directly in RStudio
boxplot(as.matrix(logCounts) ~ col(logCounts), 
        key=F, trace="none", col=Colour, margin=c(5, 5), 
        xlab="Samples", 
        ylab="Counts", main="Counts for small and large", cex.main=0.8, cex.sub=0.8, cex.lab=0.8, cex.axis=0.8)

#Visualizing the logged counts data - DensityPlot
plot(density(logCounts[,1]), ylim=c(0,1.2), xlim=c(0.2,2), col=Colour[1], margin=c(5, 5), main="Counts for small and large", cex.main=0.8, cex.sub=0.8, cex.lab=0.8, cex.axis=0.8)
for(i in 2:ncol(logCounts)) lines(density(logCounts[,i]), col=Colour[i])

#Normalizing the values in each group
countsNormalised <- scale(counts)
head(countsNormalised)

#Normalized Heatmap on samples groups
Matrix <- cor(countsNormalised)
ggcorrplot(Matrix, colors=c("goldenrod1", "white", "lightgoldenrod1"), outline.color = "white")

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Colour, label = rownames(pca_scores))) +
  geom_point(size = 2) + 
  geom_text_repel(size = 2, box.padding = 0.25, max.overlaps = 5, color = "black") +  # Make labels black
  stat_ellipse(type = "norm", linetype = 2, aes(group = Group)) +  # Add ellipses by Group
  labs(title = "Principal component analysis (PCA) ", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, color = "black", size = 0.1) +  
  geom_vline(xintercept = 0, color = "black", size = 0.1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 6, face = "bold"),  # Make axis titles bold
    axis.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 6),
    panel.grid.major = element_line(color = "black", size = 0.1),  # Major grid lines
    panel.grid.minor = element_line(color = "black", size = 0.1), # Minor grid lines
    panel.grid.major.x = element_line(color = "black", size = 0.025),
    panel.grid.minor.x = element_line(color = "black", size = 0.025),
    panel.grid.major.y = element_line(color = "black", size = 0.025),
    panel.grid.minor.y = element_line(color = "black", size = 0.025),
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +  # Create x-axis ticks
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +  # Create y-axis ticks
  scale_color_identity()  # Use the exact colors from the Colour vector

# Save the plot with high resolution
ggsave('PCA.tif', units="in", width=2.5, height=2.1, dpi=1200, compression = 'none')
