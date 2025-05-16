#Load necessary libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(svglite)
library(readxl)

#Read and process synonymous SNPs
df <- read.csv("merged_k2rt_blastp_uniprot.csv") %>%
  select(c(5:10, 21, 25, 28)) %>%
  select(Entry.Name, ssciname, Length, SNP_position, Codon_1, Codon_2, Amino_acid_1, Amino_acid_2, Is_not_synonymous) %>%
  separate(Entry.Name, into = c("Gene", "Organism"), sep = "_") %>%
  select(-Organism) %>%
  as_tibble() %>%
  select(-Is_not_synonymous) %>%
  mutate(transition = paste0(Codon_1, " (", Amino_acid_1, ") ", " → ", Codon_2, " (", Amino_acid_2, ") "),
         n.count = 1) %>%
  select(Gene, ssciname, transition, n.count) %>%
  group_by(Gene, ssciname, transition) %>%
  summarise(n.sum = sum(n.count), .groups = "drop") %>%
  pivot_wider(names_from = transition, values_from = n.sum, values_fill = list(n.sum = 0))

#Calculate column sums and filter out those with sums less than 2
column_sums <- colSums(df[ , sapply(df, is.numeric)])
columns_to_keep <- names(column_sums[column_sums >= 2])
df <- df %>%
  select(c("Gene", "ssciname", all_of(columns_to_keep)))

#Sum rows and filter those with sum greater than 2
df <- df %>%
  rowwise() %>%
  mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  filter(row_sum >= 2) %>%
  ungroup() %>%
  select(-row_sum)

#Save ssciname for splitting before converting to matrix
sscinames <- df$ssciname

#Add ID and arrange by ssciname
df <- df %>%
  mutate(ID = paste(Gene, ssciname, sep = "_")) %>%
  arrange(ssciname)

#Convert to matrix
df.matrix <- df %>%
  select(-c("Gene", "ssciname")) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()

split.s <- factor(rownames(df.matrix), levels = rownames(df.matrix))

#heat map-square
#Load the Excel sheet into R
df.ns <- read_excel("non-synonymous.xlsx")

#Create a matrix for the heatmap and remove non-numeric columns, but keep 'ID' for row names
df.ns.temp1 <- df.ns %>%
  select(-Gene, -ssciname) %>%
  column_to_rownames("ID") %>%
  as.matrix()

#Extract gene names and species names to use for row labels and grouping
gene_names <- df.ns$Gene[match(rownames(df.ns.temp1), df.ns$ID)]
species <- df.ns$ssciname[match(rownames(df.ns.temp1), df.ns$ID)]

#Clean species names (e.g., remove extra spaces and convert to lowercase for consistency)
species <- factor(str_trim(species))
genes <- factor(str_trim(gene_names))

#Color coding for the heatmap
col.ns <- colorRamp2(
  c(0, max(df.ns.temp1, na.rm = TRUE) / 5, max(df.ns.temp1, na.rm = TRUE) / 3, 
    max(df.ns.temp1, na.rm = TRUE) / 2, 2 * max(df.ns.temp1, na.rm = TRUE) / 3, max(df.ns.temp1, na.rm = TRUE)),
  c("grey99", "grey80", "grey60", "darkgoldenrod3", "goldenrod1", "lightgoldenrod1")
)

#Create the heatmap and save it as a TIFF file
svg("non-syn.svg", width = 17, height = 25, bg = "transparent", family = "arial")
Heatmap(df.ns.temp1, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = col.ns,
        row_split = species,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_labels = genes,
        row_names_side = "left",  # Move row names to the left
        row_names_rot = 0,
        heatmap_legend_param = list(title = "Counts"),
        rect_gp = gpar(col = "grey", lwd = 0.5),
        border_gp = gpar(col = "black", lty = 0.5),
        row_gap = unit(3, "mm"),
        row_title_rot = 0)
dev.off()


#heat map-square
#Load the Excel sheet into R
df.s <- read_excel("synonymous.xlsx")

#Add ID and arrange by ssciname
df.s <- df.s %>%
  mutate(ID = paste(Gene, ssciname, sep = ":")) %>%
  arrange(ssciname)

#Create a matrix for the heatmap and remove non-numeric columns, but keep 'ID' for row names
df.s.temp1 <- df.s %>%
  select(-Gene, -ssciname) %>%
  column_to_rownames("ID") %>%
  as.matrix()

#Extract gene names and species names to use for row labels and grouping
gene_names <- df.s$Gene[match(rownames(df.s.temp1), df.s$ID)]
species <- df.s$ssciname[match(rownames(df.s.temp1), df.s$ID)]

#Clean species names (e.g., remove extra spaces and convert to lowercase for consistency)
species <- factor(str_trim(species))
genes <- factor(str_trim(gene_names))

#Color coding for the heatmap
col.s <- colorRamp2(
  c(0, max(df.s.temp1, na.rm = TRUE) / 5, max(df.s.temp1, na.rm = TRUE) / 3, 
    max(df.s.temp1, na.rm = TRUE) / 2, 2 * max(df.s.temp1, na.rm = TRUE) / 3, max(df.s.temp1, na.rm = TRUE)),
  c("grey99", "grey80", "grey60", "darkgoldenrod3", "goldenrod1", "lightgoldenrod1")
)

#Create the heatmap and save it as a TIFF file
svg("syn.svg", width = 9, height = 21, bg = "transparent", family = "arial")
Heatmap(df.s.temp1, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = col.ns,
        row_split = species,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_labels = genes,
        row_names_side = "right",  # Row names on the right
        row_names_rot = 0,
        row_title_side = "right",  # Move species name (row split title) to the right
        heatmap_legend_param = list(title = "Counts"),
        rect_gp = gpar(col = "grey", lwd = 0.5),
        border_gp = gpar(col = "black", lty = 0.5),
        row_gap = unit(3, "mm"),
        row_title_rot = 0)
dev.off()