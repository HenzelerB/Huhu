#Loading the libraries for the run
library(tidyverse)
library(ggVennDiagram)
library(gplots)

#Loading the data needed for the R
COG <- read.table('COG_best.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")
KEGG <- read.table('KEGG_Annotation.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")
NR <- read.table('NR_Annotation.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")
SWISSPROT <- read.table('Swissprot_Annotation.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")

#Filter and keep only the Gene_ID column
COG_venn<- select(COG, Gene_ID)
KEGG_venn <- select(KEGG, Gene_ID)
NR_venn <- select(NR, Gene_ID)
SWISSPROT_venn <- select(SWISSPROT, Gene_ID)
gene_lists <- list(COG = COG_venn$Gene_ID, KEGG = KEGG_venn$Gene_ID, NR = NR_venn$Gene_ID, SWISSPROT = SWISSPROT_venn$Gene_ID)

#venn-diagram for the annotation
pdf('Annotation_venn.pdf', width=13, height=9, bg = "transparent")
ggVennDiagram(gene_lists, set_color = c("black", "black", "black", "black")) + scale_fill_gradient(low = "ivory", high = "goldenrod1")
dev.off()

#Loading the data needed for the R
GO <- read.table('GO_Annotation.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")
INTERPRO <- read.table('INTERPRO_Annotation.txt', sep = '\t', header = TRUE, fill = TRUE, quote = "")

#Filter and keep only the Gene_ID column
GO_venn<- select(GO, Gene_ID)
INTERPRO_venn <- select(INTERPRO, Protein_Accession)
protein_lists <- list(GO = GO_venn$Gene_ID, INTERPRO = INTERPRO_venn$Protein_Accession)

#venn-diagram for the annotation
pdf('Functional_venn.pdf', width=5, height=5, bg = "transparent")
ggVennDiagram(protein_lists, set_color = c("black", "black")) + scale_fill_gradient(low = "ivory", high = "goldenrod1")
dev.off()

#analyse and group the datatset
df.interpro <- INTERPRO %>%
  select(`Signature_Description`, `InterPro_Annotations`, `Analysis`) %>%
  group_by(`Analysis`)%>%
  summarise(Count = n())
write.csv(df.interpro, "interpro-groups.csv")