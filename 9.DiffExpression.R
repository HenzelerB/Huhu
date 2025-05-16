#Loading libraries for the run
library(tidyverse)
library(DESeq2)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(EnhancedVolcano)
library(ggplot2)
library(svglite)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(grid)
library(beeswarm)
library(VennDiagram)
library(ggVennDiagram)
library(extrafont)

#Function to abbreviate species names
abbreviate_species <- function(description) 
{species_name <- gsub(".*\\[(.*)\\].*", "\\1", description)
  abbreviated_species <- gsub("(\\w)\\w*\\s(\\w+)", "\\1.\\2", species_name)
  description <- gsub(species_name, abbreviated_species, description)
  return(description)}

#Read the data file
df <- read.csv("all.fpkm_anno.xls.csv")
df$NR_Description <- sapply(df$NR_Description, abbreviate_species)

#Summarize the data, apply ceiling, and merge back the description
df.tmp <- df %>%
  dplyr::select(NR_protein_ID,
                NR_Description,
                Huhu.S1.001_expected_count, 
                Huhu.S1.002_expected_count, 
                Huhu.S1.003_expected_count, 
                Huhu.L1.001_expected_count, 
                Huhu.L1.002_expected_count, 
                Huhu.L1.003_expected_count) %>%
  dplyr::rename(Huhu.S1 = Huhu.S1.001_expected_count,
                Huhu.S2 = Huhu.S1.002_expected_count,
                Huhu.S3 = Huhu.S1.003_expected_count,
                Huhu.L1 = Huhu.L1.001_expected_count,
                Huhu.L2 = Huhu.L1.002_expected_count,
                Huhu.L3 = Huhu.L1.003_expected_count) %>%
  group_by(NR_Description) %>%
  summarise(across(starts_with("Huhu"), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(across(starts_with("Huhu.S"), ceiling)) %>%
  mutate(across(starts_with("Huhu.L"), ceiling)) %>%
  column_to_rownames(var = "NR_Description")

#deseq2
#Differential expression setup
conds = c("Huhu.S","Huhu.S","Huhu.S","Huhu.L","Huhu.L","Huhu.L")
dds <- DESeqDataSetFromMatrix(countData = as.matrix(df.tmp), 
                              colData = data.frame(conds = factor(conds)), 
                              design = formula(~conds))
dds.tmp <- counts(dds)
write.csv(dds.tmp,"counts_deseq2.csv")

#Fit DESeq model to identify DE transcripts
dds = DESeq(dds)
colData(dds)
res = DESeq2::results(dds)
head(res)
knitr::kable(res[1:6,])
summary(res) 
res = na.omit(res)
resPadj=res[res$padj <= 0.05 , ]
dim(resPadj)
sum(res$padj <= 0.05)
sum(p.adjust(res$pvalue, method="fdr") <= 0.05)
sum(p.adjust(res$pvalue, method="holm") <= 0.01)
res = res[order(res$padj),]

#export the data outside
write.csv(as.data.frame(res),file='deseq2.csv')

#Regularizing log transformation for clustering/heat maps, etc
rld=rlogTransformation(dds)
head(assay(rld))

#Histogram to plot the relative transformed expression
hist(assay(rld), xlab="Relative distribution", ylab="Frequency",col="goldenrod2", ylim=c(0,1.5e+05), xlim=c(-5,15), main="Large vs Small", cex.main=0.8, cex.sub=0.8, cex.lab=0.8, cex.axis=0.8)

#Sample distance plot for difference between group
mycols <- brewer.pal(8, "Dark2")[1:length(unique(conds))]
sampleDist <- as.matrix(dist(t(assay(rld)), method = "euclidean")) * 0.25
svg('Sample_distance.svg', width=6, height=6)
par(family="Arial")
heatmap.2(as.matrix(sampleDist), key=FALSE, trace="none", density.info="histogram",
          col=colorpanel(200, "goldenrod2", "ivory"),
          ColSideColors=mycols[conds], RowSideColors=mycols[conds],
          margin=c(8, 8), main="Sample distance plot",
          cex.main=2, cex.sub=1.5, cex.lab=1.5, cex.axis=1.5, 
          cexRow=1.5, cexCol=1.5)
dev.off()

#Significant - Top50
top_genes_deseq2 <- as.data.frame(res) %>%
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::arrange(pvalue, desc(abs(log2FoldChange))) %>%
  head(50)
write.csv(top_genes_deseq2,"top_deseq2.csv")

#filter for A.glabripennis
res_deseq2 <- read.csv("deseq2.csv")
top_A.glabripennis <-as.data.frame(res_deseq2) %>%
  filter(grepl("\\[A\\.glabripennis\\]$", ID)) %>%
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange) > 2) %>%
  head(50)
write.csv(top_A.glabripennis,"top_A.glabripennis.csv")

#import the analysed data again
res <- read.csv("Updated_DESeq2_Data_with_Labels.csv") %>%
  select(-X)

#VolcanoPlot-DE_Transcripts
svg('VPlot_deseq2.svg', width=8, height=8)
EnhancedVolcano(res,
                lab = res$Label,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = 'p-value',
                boxedLabels = TRUE,
                xlim = c(-10, 10),
                title = 'Small vs Large',
                subtitle = "Differential expression - DESeq2",
                pCutoff = 10e-12,
                FCcutoff = 2,
                pointSize = 2,
                labSize = 4,
                col=c('grey', 'grey', 'grey', 'goldenrod1'),
                colAlpha = 0.75,
                cutoffLineType = 'dashed',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendLabels=c('Not-sig','Log2FC','p-value', 'p-value & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'black',
                max.overlaps = 30,
                axisLabSize = 12,
                titleLabSize = 14)
dev.off()

#Heat map-DE_Transcripts
top_A_glabripennis <- read.csv("top_A.glabripennis.csv")
counts_deseq2 <- read.csv("counts_deseq2.csv")
filtered_counts <- counts_deseq2 %>%
  filter(ID %in% top_A_glabripennis$ID)%>%
  mutate(ID = gsub("\\[A\\.glabripennis\\]", "", ID))
write.csv(filtered_counts, "filtered_counts_deseq2.csv", row.names = FALSE)

#DE-heatmap
hmap.DE <- read.csv("hmap-edited.csv")%>%
  mutate(across(-ID, ~ (. - min(.)) / (max(.) - min(.))))%>%
  column_to_rownames("ID") %>%
  as.matrix()
#Color coding for the heatmap
col.s <- colorRamp2(c(0, max(hmap.DE, na.rm = TRUE) / 4, max(hmap.DE, na.rm = TRUE) / 3, max(hmap.DE, na.rm = TRUE) / 2, max(hmap.DE, na.rm = TRUE)), c("ivory", "lightgoldenrod2", "springgreen3", "palegreen3", "salmon"))
svg("DE_Hmap.svg", width = 4, height = 8, bg = "transparent", family = "arial")
Heatmap(hmap.DE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col.s,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left",
        column_names_rot = 90,  
        column_title = "A.glabripennis normalised counts",  
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(title = "Counts"),
        rect_gp = gpar(col = "grey", lwd = 0.5),
        border_gp = gpar(col = "black", lty = 0.5),
        column_gap = unit(2, "mm"),
        column_title_rot = 0,
        # Add cell values as text within each heatmap block
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", hmap.DE[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        })
dev.off()

#for PPI
PPI <- df %>%
  filter(grepl("\\[Anoplophora glabripennis\\]$", NR_Description))

#limma
#read dataset for limma and do the needed analysis
df.tmp <- df %>%
  dplyr::select(NR_protein_ID,
                NR_Description,
                Huhu.S1.001_expected_count, 
                Huhu.S1.002_expected_count, 
                Huhu.S1.003_expected_count, 
                Huhu.L1.001_expected_count, 
                Huhu.L1.002_expected_count, 
                Huhu.L1.003_expected_count) %>%
  dplyr::rename(Huhu.S1 = Huhu.S1.001_expected_count,
                Huhu.S2 = Huhu.S1.002_expected_count,
                Huhu.S3 = Huhu.S1.003_expected_count,
                Huhu.L1 = Huhu.L1.001_expected_count,
                Huhu.L2 = Huhu.L1.002_expected_count,
                Huhu.L3 = Huhu.L1.003_expected_count) %>%
  group_by(NR_Description) %>%
  summarise(across(starts_with("Huhu"), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(across(starts_with("Huhu.S"), ceiling)) %>%
  mutate(across(starts_with("Huhu.L"), ceiling))
counts <- as.data.frame(df.tmp[, 2:7])
rownames(counts) <- df.tmp$NR_Description
head(counts)
colSums(counts)
colSums(counts == 0)

#process the dge list
dge = DGEList(counts=counts)
dge = calcNormFactors(dge)
logCPM = cpm(dge, log=TRUE, prior.count=3)
options(width=100)
head(logCPM, 3)

#The beeswarm
t.test(logCPM[6,] ~ conds)
design = model.matrix(~conds)
design
v = voom(dge, design, plot=TRUE)
q = voom(dge, design, plot=TRUE, normalize="quantile")

#statistical testing
fit = lmFit(v, design)
fit = eBayes(fit)
tt = topTable(fit, coef=ncol(design), n=nrow(counts))
head(tt)
sum(tt$adj.P.Val < 0.05)
sum(p.adjust(tt$P.Value, method="fdr") <= 0.05)

#export the data outside
write.csv(as.data.frame(tt),file='limma.csv')

#volcano plot
volcanoplot(fit, coef=2)
sigGenes = which(tt$adj.P.Val <= 0.05)
length(sigGenes)
volcanoplot(fit, coef=2)
points(tt$logFC[sigGenes], -log10(tt$P.Value[sigGenes]), col='red', pch=16, cex=0.5)
sigGenes = which(tt$adj.P.Val <= 0.05 & (abs(tt$logFC) > log2(2)))
length(sigGenes)
volcanoplot(fit, coef=2)
points(tt$logFC[sigGenes], -log10(tt$P.Value[sigGenes]), col='red', pch=16, cex=0.5)
limmaPadj = tt[tt$adj.P.Val <= 0.05, ]

#edgeR
#read dataset for edgeR and do the needed analysis
df.tmp <- df %>%
  dplyr::select(NR_protein_ID,
                NR_Description,
                Huhu.S1.001_expected_count, 
                Huhu.S1.002_expected_count, 
                Huhu.S1.003_expected_count, 
                Huhu.L1.001_expected_count, 
                Huhu.L1.002_expected_count, 
                Huhu.L1.003_expected_count) %>%
  dplyr::rename(Huhu.S1 = Huhu.S1.001_expected_count,
                Huhu.S2 = Huhu.S1.002_expected_count,
                Huhu.S3 = Huhu.S1.003_expected_count,
                Huhu.L1 = Huhu.L1.001_expected_count,
                Huhu.L2 = Huhu.L1.002_expected_count,
                Huhu.L3 = Huhu.L1.003_expected_count) %>%
  group_by(NR_Description) %>%
  summarise(across(starts_with("Huhu"), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  mutate(across(starts_with("Huhu.S"), ceiling)) %>%
  mutate(across(starts_with("Huhu.L"), ceiling))
counts <- as.data.frame(df.tmp[, 2:7])
rownames(counts) <- df.tmp$NR_Description
head(counts)
colSums(counts)
colSums(counts == 0)

#DGEList for the edge run
y = DGEList(counts=counts, group=conds)
y = calcNormFactors(y)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)
et = exactTest(y) 
knitr::kable(topTags(et, n=4)$table)
edge = as.data.frame(topTags(et, n=nrow(counts))) 
sum(edge$FDR <= 0.05)
sum(p.adjust(edge$PValue, method="fdr") <= 0.05)
edgePadj = edge[edge$FDR <= 0.05, ]

#export the data outside
write.csv(as.data.frame(edge),file='edge.csv')

#Common genes between the 3 different packages
#Venn diagram 
setlist = list(edgeR=rownames(edgePadj), DESeq2=rownames(resPadj))
venn(setlist)
setlist = list(edgeR=rownames(edgePadj), 
               DESeq2=rownames(resPadj),
               Limma=rownames(limmaPadj))
pdf("DE_venn.pdf", width = 5, height = 5, bg = "transparent")
ggVennDiagram(setlist, set_color = c("black", "black", "black")) + scale_fill_gradient(low = "ivory", high = "goldenrod1")
dev.off()

#merge and split the dataframe to create grouped dataset for the remaining analysis
deseq2 <- read.csv("deseq2.csv")
edgeR <- read.csv("edge.csv")
limma <- read.csv("limma.csv")

#merge the datasets based on common IDs
df.m <- deseq2 %>%
  inner_join(edgeR %>% select(ID), by = "ID") %>%
  inner_join(limma %>% select(ID), by = "ID") %>%
  select(names(deseq2)) %>%
  mutate(Species = sub(".*\\[(.*)\\].*", "\\1", ID))
write.csv(df.m, "merged.csv", row.names = FALSE)

#Significant_deseq
top_genes_deseq2.m <- as.data.frame(df.m) %>%
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::arrange(pvalue, desc(abs(log2FoldChange))) %>%
  head(100)
write.csv(top_genes_deseq2.m,"top_merged.csv")

#filter for A.glabripennis
res_deseq2.m <- read.csv("merged.csv") %>%
  as.data.frame()
top_A.glabripennis.m <-as.data.frame(res_deseq2.m) %>%
  filter(grepl("\\[A\\.glabripennis\\]$", ID)) %>%
  dplyr::filter(pvalue < 0.05 & abs(log2FoldChange) > 2) %>%
  head(100)
write.csv(top_A.glabripennis.m,"top_A.glabripennis_merged.csv")

#import the analysed data again adn filter it
res.M <- read.csv("top_deseq2_M_labelled.csv")
merged_with_labels <- res_deseq2.m %>%
  left_join(res.M %>% select(ID, Label), by = "ID") %>%
  mutate(Gene.ID = gsub("\\[.*\\]", "", Label) %>% trimws())
write.csv(merged_with_labels,"merged_DE_labelled.csv")

#VolcanoPlot-DE_Transcripts
svg('VPlot_deseq2.svg', width=8, height=9)
EnhancedVolcano(merged_with_labels,
                lab = merged_with_labels$Gene.ID,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = 'p-value',
                boxedLabels = TRUE,
                xlim = c(-10, 10),
                title = 'Small vs Large',
                subtitle = "Differential expression - DESeq2",
                pCutoff = 10e-12,
                FCcutoff = 2,
                pointSize = 1,
                labSize = 3,
                col=c('grey', 'grey', 'grey', 'goldenrod1'),
                colAlpha = 0.75,
                cutoffLineType = 'dashed',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendLabels=c('Not-sig','Log2FC','p-value', 'p-value & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 12,
                drawConnectors = TRUE,
                widthConnectors = 0.15,
                colConnectors = 'black',
                max.overlaps = 75,
                axisLabSize = 12,
                titleLabSize = 14)+
  scale_x_continuous(breaks = seq(-10, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 300, by = 25))
dev.off()


