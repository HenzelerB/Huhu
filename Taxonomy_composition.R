# =============================================================================
# Taxonomic composition of the DE set + assembly-level contamination screen
# -----------------------------------------------------------------------------
# Purpose (addresses Reviewer 4, comment 5):
#   Quantify how non-arthropod / unassigned transcripts affect the results.
#   Produces:
#     DE_taxonomy_composition.csv     - every DE transcript + best-hit species + class
#     DE_nonarthropod_transcripts.csv - just the non-arthropod DE transcripts
#     Figure_Contamination.pdf/.png   - (A) assembly composition (by transcript vs
#                                        coverage, from BlobTools); (B) DE-set
#                                        composition; (C) non-arthropod breakdown
#
# IMPORTANT — why classification is by best-hit species, not by joining BlobTools:
#   The BlobTools assembly and the DE/quantification assembly are DIFFERENT Trinity
#   runs (Trinity numbers contigs non-deterministically), so their DN identifiers
#   do NOT correspond (~30% string overlap only). A per-transcript ID join is
#   therefore invalid. Instead we classify each DE transcript by the best-hit
#   species already embedded in its NR_Description ("... [Genus species]"), which
#   is the annotation that actually drove functional interpretation.
#
# HOW TO RUN:
#   1. set DE_FILE to your intersection DE table (a CSV whose description column
#      contains "[Genus species]"; e.g. sensitivity_reference_DE.csv, or Table S5).
#   2. set BLOB_FILE to the BlobTools table (…blobDB.table.txt).
#   3. Rscript Taxonomy_composition.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

# ------------------------------- CONFIG --------------------------------------
DE_FILE   <- "sensitivity_out/sensitivity_reference_DE.csv"  # DE set (intersection)
BLOB_FILE <- file.choose()   # easiest: pick the blobDB.table.txt file in the dialog
OUTDIR    <- "sensitivity_out"
dir.create(OUTDIR, showWarnings = FALSE)

# Colours (match the sensitivity figure)
amber <- "#E0A100"; grey <- "#8A8A8A"; ink <- "#222222"; lgrey <- "#CFCFCF"; red <- "#B00020"

# -----------------------------------------------------------------------------
# 1. Classification helper: best-hit species -> arthropod / non-arthropod
# -----------------------------------------------------------------------------
# Arthropod genera observed in this dataset. This list is explicit and editable;
# a genus not listed here is classed "non-arthropod" and will appear in the
# non-arthropod CSV, so any additions are easy to spot and add.
# (For a fully taxonomy-driven alternative see the commented taxize block below.)
ARTHROPOD_GENERA <- c(
  "Abscondita", "Acromyrmex", "Aedes", "Aethina", "Agrilus", "Alphitobius", "Anaspis",
  "Anatolica", "Anopheles", "Anoplophora", "Aphis", "Apis", "Apolygus", "Apriona",
  "Argiope", "Asbolus", "Athalia", "Bactrocera", "Bemisia", "Biphyllus", "Blattella",
  "Bombyx", "Bradysia", "Camponotus", "Centruroides", "Ceratitis", "Ceratosolen", "Chilo",
  "Chrysomela", "Cinara", "Clunio", "Colaphellus", "Coleoptera", "Contarinia",
  "Coptotermes", "Cryptotermes", "Ctenocephalides", "Culex", "Cyclommatus", "Cyphomyrmex",
  "Danaus", "Daphnia", "Dendroctonus", "Diabrotica", "Dinoponera", "Drosophila", "Ephemera",
  "Epicauta", "Eumeta", "Fopius", "Frankliniella", "Galleria", "Gnatocerus", "Graminella",
  "Halyomorpha", "Harmonia", "Harpegnathos", "Heliothis", "Hermetia", "Hyalella", "Hycleus",
  "Hymaea", "Hypothenemus", "Ignelater", "Ischnomera", "Lamprigera", "Leptinotarsa",
  "Limulus", "Linepithema", "Lytta", "Manduca", "Megachile", "Meladema", "Melipona",
  "Meloe", "Mesosa", "Microdera", "Monochamus", "Monomorium", "Mycetophagus", "Mylabris",
  "Nacerdes", "Nasonia", "Nicrophorus", "Nilaparvata", "Nymphon", "Odontomachus",
  "Oedemera", "Oedemeridae", "Onthophagus", "Ooceraea", "Oppiella", "Orussus", "Oryctes",
  "Osmia", "Ostrinia", "Pachyrhynchus", "Papilio", "Pediculus", "Penaeus", "Perophthalma",
  "Photinus", "Plutella", "Pseudoatta", "Reticulitermes", "Rhagoletis", "Rhipicephalus",
  "Rhynchophorus", "Riptortus", "Sarcoptes", "Sitophilus", "Stegodyphus", "Stilpnonotus",
  "Teleopsis", "Temnothorax", "Tenebrio", "Tenebrionoidea", "Tetranychus", "Thrips",
  "Tigriopus", "Timarcha", "Timema", "Trachymyrmex", "Tribolium", "Trichogramma",
  "Trichomalopsis", "Trictenotoma", "Trinorchestia", "Wasmannia", "Xylotrechus", "Zerene",
  "Zootermopsis", "Zygaena")

# Sub-categories for the non-arthropod best hits (used in figure panel C)
APICOMPLEXA <- c("Gregarina","Plasmodium","Cryptosporidium","Eimeria","Toxoplasma","Neospora",
                 "Besnoitia","Theileria","Babesia","Cyclospora","Cardiosporidium","Perkinsus")
BACTERIA    <- c("Yersinia","Edwardsiella","Proteus","Morganella")
OTHERPROT   <- c("Stylonychia","Tetrahymena","Trypanosoma","Naegleria","Capsaspora",
                 "Acytostelium","Planoprotostelium","Tieghemostelium","Paulinella")
NEMATODE    <- c("Trichinella","Aphelenchoides","Hymenolepis","Diploscapter")

genus_of  <- function(sp) ifelse(is.na(sp) | sp=="", NA_character_, sub("\\s.*$", "", sp))
classify  <- function(sp) ifelse(genus_of(sp) %in% ARTHROPOD_GENERA, "arthropod", "non-arthropod")
subgroup  <- function(sp) {
  g <- genus_of(sp)
  dplyr::case_when(
    g %in% APICOMPLEXA ~ "Apicomplexan gut parasite",
    g %in% BACTERIA    ~ "Gut bacteria",
    g %in% OTHERPROT   ~ "Other protist",
    g %in% NEMATODE    ~ "Nematode/helminth",
    TRUE               ~ "Vertebrate/other metazoan (conserved domains)")
}

# ---- (optional) rigorous taxonomy-based classification via NCBI --------------
# If you prefer to resolve phylum from NCBI taxonomy instead of the genus list:
#   library(taxize)   # needs internet / an ENTREZ key
#   uniq <- unique(de$best_hit_species)
#   cl   <- taxize::classification(uniq, db = "ncbi")
#   is_arth <- sapply(cl, function(x) "Arthropoda" %in% x$name)
#   ... then map back to de$best_hit_species.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2. Load DE set, extract best-hit species, classify
# -----------------------------------------------------------------------------
de <- read.csv(DE_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# auto-detect the description column (character col with the most "[...]" tags)
desc_col <- names(de)[which.max(sapply(de, function(col)
  if (is.character(col)) sum(grepl("\\[.*\\]", col)) else -1))]
lfc_col  <- names(de)[grep("lfc|log2|fold", names(de), ignore.case = TRUE)][1]
fdr_col  <- names(de)[grep("fdr|padj|adj", names(de), ignore.case = TRUE)][1]

de <- de %>%
  mutate(NR_Description   = .data[[desc_col]],
         best_hit_species = sub(".*\\[([^]]+)\\]\\s*$", "\\1", NR_Description),
         best_hit_species = ifelse(grepl("\\[", NR_Description), best_hit_species, NA),
         class            = classify(best_hit_species))

out <- de %>% transmute(NR_Description,
                        log2FC = if (!is.na(lfc_col)) .data[[lfc_col]] else NA,
                        padj   = if (!is.na(fdr_col)) .data[[fdr_col]] else NA,
                        best_hit_species, class)
write.csv(out, file.path(OUTDIR, "DE_taxonomy_composition.csv"), row.names = FALSE)
write.csv(filter(out, class == "non-arthropod"),
          file.path(OUTDIR, "DE_nonarthropod_transcripts.csv"), row.names = FALSE)

n_de   <- nrow(out)
n_arth <- sum(out$class == "arthropod")
n_non  <- n_de - n_arth
n_unan <- sum(is.na(out$best_hit_species))
cat(sprintf("DE transcripts: %d | arthropod %d (%.1f%%) | non-arthropod %d (%.1f%%) | unassigned %d\n",
            n_de, n_arth, 100*n_arth/n_de, n_non, 100*n_non/n_de, n_unan))

# -----------------------------------------------------------------------------
# 3. Assembly-level composition from the BlobTools table
# -----------------------------------------------------------------------------
lines <- readLines(BLOB_FILE)
lines <- lines[!grepl("^##", lines)]                 # drop ## metadata
hdr   <- strsplit(sub("^#\\s*", "", lines[1]), "\t")[[1]]
blob  <- read.table(text = lines[-1], sep = "\t", quote = "", comment.char = "",
                    stringsAsFactors = FALSE)
# columns are positional: 11 = cov_sum, 12 = phylum
blob  <- data.frame(cov_sum = as.numeric(blob[[11]]), phylum = blob[[12]],
                    stringsAsFactors = FALSE)

cat_of <- function(phy) dplyr::case_when(
  phy == "Arthropoda" ~ "Arthropoda",
  phy %in% c("no-hit","undef","cellular organisms-undef") ~ "Unassigned/no-hit",
  phy == "Chordata" ~ "Chordata (conserved-domain)",
  phy %in% c("Ascomycota","Basidiomycota") ~ "Fungi",
  phy == "Apicomplexa" ~ "Protozoa (Apicomplexa)",
  grepl("bacter|Pseudomonadota|Bacteroidota|Spirochaetota|Planctomycetota|Campylobacterota|Bacillota|Actinomycetota|Myxococcota|Nitrospirota|Acidobacteriota|Chlamydiota|Cyanobacteriota|Deinococcota|Firmicutes|Proteobacteria",
        phy, ignore.case = TRUE) ~ "Bacteria",
  TRUE ~ "Other eukaryote")

cat_levels <- c("Arthropoda","Chordata (conserved-domain)","Bacteria","Fungi",
                "Protozoa (Apicomplexa)","Other eukaryote","Unassigned/no-hit")
asm <- blob %>%
  mutate(cat = factor(cat_of(phylum), levels = cat_levels)) %>%
  group_by(cat) %>%
  summarise(n = n(), cov = sum(cov_sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(`% of transcripts` = 100*n/sum(n),
         `% of coverage (reads)` = 100*cov/sum(cov))

# -----------------------------------------------------------------------------
# 4. Figure
# -----------------------------------------------------------------------------
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        axis.title = element_text(size = 9), legend.position = "bottom")

# Panel A: assembly composition (grouped bars)
pA <- asm %>%
  select(cat, `% of transcripts`, `% of coverage (reads)`) %>%
  pivot_longer(-cat, names_to = "measure", values_to = "pct") %>%
  ggplot(aes(cat, pct, fill = measure)) +
  geom_col(position = position_dodge(.75), width = .68) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(.75),
            hjust = -0.1, size = 2.7) +
  scale_fill_manual(values = c("% of transcripts" = grey, "% of coverage (reads)" = amber)) +
  coord_flip(ylim = c(0, 68)) + scale_x_discrete(limits = rev(cat_levels)) +
  labs(title = sprintf("A  Assembly composition (BlobTools, %s transcripts)",
                       format(nrow(blob), big.mark = ",")),
       x = NULL, y = "percentage of assembly", fill = NULL) + th

# Panel B: DE-set composition (stacked)
pB <- data.frame(x = "DE set",
                 class = factor(c("arthropod","non-arthropod"),
                                levels = c("non-arthropod","arthropod")),
                 pct = c(100*n_arth/n_de, 100*n_non/n_de)) %>%
  ggplot(aes(x, pct, fill = class)) +
  geom_col(width = .5) +
  geom_text(data = data.frame(x="DE set", y=100*n_arth/n_de/2,
                              lab=sprintf("arthropod\n%.1f%%", 100*n_arth/n_de)),
            aes(x, y, label = lab), inherit.aes = FALSE, colour = "white",
            fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c(arthropod = amber, `non-arthropod` = red)) +
  coord_flip() +
  labs(title = sprintf("B  DE set (n=%s): %.1f%% arthropod, %.1f%% non-arthropod, %d unassigned",
                       format(n_de, big.mark=","), 100*n_arth/n_de, 100*n_non/n_de, n_unan),
       x = NULL, y = "percent of DE set", fill = NULL) + th

# Panel C: breakdown of the non-arthropod best hits
pC <- out %>% filter(class == "non-arthropod") %>%
  mutate(grp = subgroup(best_hit_species)) %>%
  count(grp) %>%
  ggplot(aes(reorder(grp, n), n)) +
  geom_col(fill = lgrey, width = .68) +
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  coord_flip() +
  labs(title = sprintf("C  Non-arthropod best hits (n=%d)", n_non),
       x = NULL, y = "transcripts") + th

fig <- pA / (pB + pC + plot_layout(widths = c(1, 1))) + plot_layout(heights = c(1.25, 1))
ggsave(file.path(OUTDIR, "Figure_Contamination.pdf"), fig, width = 11, height = 8)
ggsave(file.path(OUTDIR, "Figure_Contamination.png"), fig, width = 11, height = 8, dpi = 200)
cat("Wrote DE_taxonomy_composition.csv, DE_nonarthropod_transcripts.csv, Figure_Contamination.pdf/.png\n")
