# =============================================================================
# Sensitivity analysis for the small-vs-large P. reticularis DE comparison
# -----------------------------------------------------------------------------
# Purpose (addresses Reviewer 4, comment 1):
#   Show that the differentially-expressed (DE) transcript set is stable and is
#   NOT driven by any single pooled replicate, by re-running the FULL DE
#   pipeline (DESeq2 AND edgeR AND limma-voom, taking the intersection, exactly
#   as in the manuscript) on:
#     (i)  leave-one-out (LOO) subsets   - drop each library in turn (6 runs)
#     (ii) 2-vs-2 balanced subsets       - every choice of 2 small + 2 large (9 runs)
#     (iii) read-depth down-sampling      - thin counts to 90..50% (optional)
#   For each subset we quantify how much of the full-data DE set is recovered
#   and how well fold-changes agree.
#
# Inputs:
#   all.fpkm_anno.xls.csv  (the same RSEM expected_count table used originally,
#   containing NR_Description + the six *_expected_count columns).
#
# Outputs (written to OUTDIR):
#   sensitivity_summary.csv          - one row per subset run, all metrics
#   sensitivity_reference_DE.csv     - the full-data reference DE set
#   Figure_Sensitivity.pdf / .png    - multi-panel supplementary figure
#   (optional) per-run DE gene lists
#
# HOW TO RUN ON YOUR DATA:
#   1. set DEMO <- FALSE
#   2. set COUNT_FILE to your all.fpkm_anno.xls.csv path
#   3. Rscript Sensitivity_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma)
  library(dplyr);  library(tidyr); library(ggplot2); library(patchwork)
})

# ------------------------------- CONFIG --------------------------------------
DEMO        <- FALSE           # TRUE = simulate data to test the script; FALSE = use your file
COUNT_FILE  <- "all_fpkm_anno_xls.csv"
OUTDIR      <- "sensitivity_out"
PADJ_CUT    <- 0.01            # manuscript threshold
LFC_CUT     <- 2              # |log2FC| threshold
TOPN        <- 50             # "top DE transcripts" set size for recovery metric
DO_DOWNSAMPLE <- TRUE          # read-depth thinning module (set FALSE to skip / speed up)
DS_DEPTHS   <- c(0.9,0.75,0.6,0.5)
DS_REPS     <- 3             # bootstrap thinnings per depth
SEED        <- 1
set.seed(SEED)
dir.create(OUTDIR, showWarnings = FALSE)

# Column mapping: original *_expected_count columns -> tidy sample IDs
# (file uses hyphens: Huhu-S1-001_expected_count etc.)
COLMAP <- c(
  Huhu.S1 = "Huhu-S1-001_expected_count",
  Huhu.S2 = "Huhu-S1-002_expected_count",
  Huhu.S3 = "Huhu-S1-003_expected_count",
  Huhu.L1 = "Huhu-L1-001_expected_count",
  Huhu.L2 = "Huhu-L1-002_expected_count",
  Huhu.L3 = "Huhu-L1-003_expected_count"
)
SMALL <- c("Huhu.S1","Huhu.S2","Huhu.S3")
LARGE <- c("Huhu.L1","Huhu.L2","Huhu.L3")

# --------------------------- DATA PREPARATION --------------------------------
abbreviate_species <- function(description) {
  species_name <- gsub(".*\\[(.*)\\].*", "\\1", description)
  abbreviated  <- gsub("(\\w)\\w*\\s(\\w+)", "\\1.\\2", species_name)
  gsub(species_name, abbreviated, description, fixed = TRUE)
}

# Reproduces the original preprocessing: collapse transcripts by NR_Description
# (mean of member transcripts), ceiling to integers. Done ONCE on all 6 columns
# so gene definitions are identical across every subset (keeps DE sets comparable).
prep_counts <- function(file) {
  # read only the needed columns for speed/memory (respecting CSV quoting)
  hdr <- names(read.csv(file, nrows = 1, check.names = FALSE))
  need <- c("NR_Description", unname(COLMAP))
  stopifnot(all(need %in% hdr))
  cc <- ifelse(hdr %in% need,
               ifelse(hdr == "NR_Description", "character", "numeric"), "NULL")
  df <- read.csv(file, check.names = FALSE, colClasses = cc)
  m <- df[, c("NR_Description", unname(COLMAP))]
  colnames(m) <- c("NR_Description", names(COLMAP))
  m %>%
    group_by(NR_Description) %>%
    summarise(across(all_of(names(COLMAP)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(all_of(names(COLMAP)), ceiling)) %>%
    as.data.frame() -> agg
  rn <- agg$NR_Description
  mat <- as.matrix(agg[, names(COLMAP)])
  rownames(mat) <- rn
  storage.mode(mat) <- "integer"
  # drop the unannotated supergroup and all-zero-count rows. No further low-count
  # pre-filter: this reproduces the original pipeline (DESeq2's own independent
  # filtering handles low-count genes), matching the published DESeq2 DEG count.
  mat <- mat[rownames(mat) != "--", , drop = FALSE]
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat
}

# Simulated matrix for DEMO/testing: NB counts, a defined subset truly DE.
simulate_counts <- function(n_genes = 4000, n_de = 400) {
  base <- rnbinom(n_genes, mu = 300, size = 2) + 5
  mk <- function(lfc) {
    sapply(base * 2^lfc, function(mu) rnbinom(1, mu = max(mu, 1), size = 3))
  }
  lfc <- rep(0, n_genes)
  de_idx <- sample(n_genes, n_de)
  lfc[de_idx] <- sample(c(-1,1), n_de, TRUE) * runif(n_de, 2.2, 5)  # large-group effect
  cols <- list()
  for (s in SMALL) cols[[s]] <- mk(rep(0, n_genes))
  for (s in LARGE) cols[[s]] <- mk(lfc)
  mat <- do.call(cbind, cols)
  rownames(mat) <- sprintf("gene_%04d", seq_len(n_genes))
  storage.mode(mat) <- "integer"
  mat
}

# ------------------------- DE ENGINES (3 tools) ------------------------------
group_of <- function(samples)
  factor(ifelse(samples %in% SMALL, "S", "L"), levels = c("S","L"))

de_deseq2 <- function(mat, samples) {
  grp <- group_of(samples)
  dds <- DESeqDataSetFromMatrix(mat[, samples, drop = FALSE],
                                data.frame(grp = grp), design = ~ grp)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("grp","L","S")))
  out <- data.frame(gene = rownames(res), lfc = res$log2FoldChange, fdr = res$padj,
                    row.names = NULL)
  rm(dds, res); gc(FALSE); out
}

de_edger <- function(mat, samples) {
  grp <- group_of(samples)
  y <- DGEList(counts = mat[, samples, drop = FALSE], group = grp)
  y <- calcNormFactors(y)
  des <- model.matrix(~ grp)
  y <- estimateDisp(y, des)
  fit <- glmQLFit(y, des)
  qlf <- glmQLFTest(fit, coef = 2)          # coef 2 = L vs S
  tt <- topTags(qlf, n = Inf)$table
  out <- data.frame(gene = rownames(tt), lfc = tt$logFC, fdr = tt$FDR, row.names = NULL)
  rm(y, fit, qlf, tt); gc(FALSE); out
}

de_limma <- function(mat, samples) {
  grp <- group_of(samples)
  y <- DGEList(counts = mat[, samples, drop = FALSE], group = grp)
  y <- calcNormFactors(y)
  des <- model.matrix(~ grp)
  v <- voom(y, des)
  fit <- eBayes(lmFit(v, des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  out <- data.frame(gene = rownames(tt), lfc = tt$logFC, fdr = tt$adj.P.Val, row.names = NULL)
  rm(y, v, fit, tt); gc(FALSE); out
}

sig_set <- function(de) de$gene[!is.na(de$fdr) & de$fdr < PADJ_CUT & abs(de$lfc) > LFC_CUT]

# Run all three tools on a sample subset; return per-tool sig sets, the
# intersection (manuscript pipeline), and the DESeq2 table (fold-change ref).
run_pipeline <- function(mat, samples) {
  d <- de_deseq2(mat, samples)
  e <- de_edger(mat, samples)
  l <- de_limma(mat, samples)
  inter <- Reduce(intersect, list(sig_set(d), sig_set(e), sig_set(l)))
  list(intersection = inter, deseq2 = d,
       n_deseq2 = length(sig_set(d)), n_edger = length(sig_set(e)),
       n_limma = length(sig_set(l)))
}

# ------------------------------- METRICS -------------------------------------
jaccard  <- function(a, b) if (length(union(a,b))==0) NA else length(intersect(a,b))/length(union(a,b))
recovery <- function(ref, sub) if (length(ref)==0) NA else length(intersect(ref,sub))/length(ref)

metrics_vs_ref <- function(label, type, res, ref_set, ref_top, ref_lfc) {
  sub <- res$intersection
  shared <- intersect(ref_lfc$gene, res$deseq2$gene)
  lfc_r <- ref_lfc$lfc[match(shared, ref_lfc$gene)]
  lfc_s <- res$deseq2$lfc[match(shared, res$deseq2$gene)]
  ok <- is.finite(lfc_r) & is.finite(lfc_s)
  data.frame(
    run = label, type = type,
    n_DEG = length(sub),
    overlap_with_ref = length(intersect(ref_set, sub)),
    recovery = recovery(ref_set, sub),           # fraction of full DE set recovered
    jaccard  = jaccard(ref_set, sub),
    top_recovery = recovery(ref_top, sub),       # fraction of full top-N recovered
    lfc_pearson  = suppressWarnings(cor(lfc_r[ok], lfc_s[ok], method = "pearson")),
    lfc_spearman = suppressWarnings(cor(lfc_r[ok], lfc_s[ok], method = "spearman")),
    n_deseq2 = res$n_deseq2, n_edger = res$n_edger, n_limma = res$n_limma
  )
}

# ------------------------------- RUN -----------------------------------------
message(">> Preparing counts ...")
counts <- if (DEMO) simulate_counts() else prep_counts(COUNT_FILE)
message(sprintf("   %d genes x %d samples", nrow(counts), ncol(counts)))

message(">> Reference (full 3 vs 3) ...")
ref <- run_pipeline(counts, c(SMALL, LARGE))
ref_set <- ref$intersection
# top-N reference transcripts by DESeq2 significance among the intersection set
ref_d <- ref$deseq2
ref_top <- ref_d %>% filter(gene %in% ref_set) %>% arrange(fdr) %>% slice_head(n = TOPN) %>% pull(gene)
message(sprintf("   reference intersection DE set: %d transcripts", length(ref_set)))
write.csv(ref_d %>% filter(gene %in% ref_set) %>% arrange(fdr),
          file.path(OUTDIR, "sensitivity_reference_DE.csv"), row.names = FALSE)

results <- list()

message(">> Leave-one-out (6 runs) ...")
for (drop in c(SMALL, LARGE)) {
  keep <- setdiff(c(SMALL, LARGE), drop)
  r <- run_pipeline(counts, keep)
  results[[paste0("LOO_drop_", drop)]] <-
    metrics_vs_ref(paste0("drop ", drop), "LOO", r, ref_set, ref_top, ref_d)
  message(sprintf("   dropped %s: %d DEG, recovery %.2f", drop,
                  length(r$intersection), recovery(ref_set, r$intersection)))
  rm(r); gc(FALSE)
}

message(">> Balanced 2-vs-2 subsets (9 runs) ...")
for (si in combn(SMALL, 2, simplify = FALSE)) {
  for (li in combn(LARGE, 2, simplify = FALSE)) {
    lab <- paste(c(si, li), collapse = "+")
    r <- run_pipeline(counts, c(si, li))
    results[[paste0("SUB_", lab)]] <-
      metrics_vs_ref(lab, "subset_2v2", r, ref_set, ref_top, ref_d)
    message(sprintf("   %s: %d DEG, recovery %.2f", lab,
                    length(r$intersection), recovery(ref_set, r$intersection)))
    rm(r); gc(FALSE)
  }
}

summary_df <- bind_rows(results)
write.csv(summary_df, file.path(OUTDIR, "sensitivity_summary.csv"), row.names = FALSE)

# ---------------------- optional read-depth downsampling ---------------------
ds_df <- NULL
if (DO_DOWNSAMPLE) {
  message(">> Read-depth down-sampling (DESeq2 recovery curve) ...")
  thin <- function(m, p) {
    matrix(rbinom(length(m), size = as.vector(m), prob = p),
           nrow = nrow(m), dimnames = dimnames(m))
  }
  ds_rows <- list()
  for (p in DS_DEPTHS) {
    for (rep in seq_len(DS_REPS)) {
      tm <- thin(counts, p); storage.mode(tm) <- "integer"
      d <- de_deseq2(tm, c(SMALL, LARGE))
      rec <- recovery(ref_set, sig_set(d))
      ds_rows[[paste(p, rep)]] <- data.frame(depth = p, rep = rep, recovery = rec)
      rm(tm, d); gc(FALSE)
    }
    message(sprintf("   depth %.0f%%: mean recovery %.2f", p*100,
                    mean(sapply(ds_rows[grepl(paste0("^",p," "), names(ds_rows))],
                                function(x) x$recovery))))
  }
  ds_df <- bind_rows(ds_rows)
  write.csv(ds_df, file.path(OUTDIR, "downsampling_recovery.csv"), row.names = FALSE)
}

# ------------------------------- FIGURE --------------------------------------
amber <- "#E0A100"; grey <- "#8A8A8A"; ink <- "#222222"; blue <- "#5A7DA0"
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9, colour = grey),
        axis.title = element_text(size = 10),
        legend.position = "bottom")
lab_type <- c(LOO = "Leave-one-out (n=3->2v3)", subset_2v2 = "2-vs-2 subset (n=2)")
summary_df$type_lab <- lab_type[summary_df$type]

# A: fold-change concordance across ALL subsets (the stability headline)
pA <- summary_df %>%
  select(type_lab, lfc_pearson, lfc_spearman) %>%
  pivot_longer(c(lfc_pearson, lfc_spearman), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, lfc_pearson = "Pearson", lfc_spearman = "Spearman")) %>%
  ggplot(aes(metric, value, fill = type_lab)) +
  geom_boxplot(outlier.shape = NA, alpha = .35, colour = grey, position = position_dodge(.8)) +
  geom_point(aes(colour = type_lab),
             position = position_jitterdodge(jitter.width = .12, dodge.width = .8),
             size = 1.7, alpha = .85) +
  scale_fill_manual(values = c(amber, grey)) + scale_colour_manual(values = c(amber, grey)) +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(title = "A  Fold-change concordance with full data",
       subtitle = "effect estimates stable in every subset (r >= 0.90)",
       x = NULL, y = "correlation of log2FC", fill = NULL, colour = NULL) + th

# B: leave-one-out recovery (the direct n=3 stability test)
pB <- summary_df %>% filter(type == "LOO") %>%
  select(run, recovery, top_recovery) %>%
  mutate(run = gsub("Huhu.|LOO_drop_", "", run)) %>%
  pivot_longer(c(recovery, top_recovery), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, recovery = "DE-set recovery", top_recovery = "Top-50 recovery")) %>%
  ggplot(aes(run, value, fill = metric)) +
  geom_col(position = position_dodge(.75), width = .7, alpha = .9) +
  scale_fill_manual(values = c("DE-set recovery" = amber, "Top-50 recovery" = blue)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "B  Leave-one-out recovery",
       subtitle = sprintf("drop each pool in turn; full set = %s DE transcripts", format(length(ref_set), big.mark=",")),
       x = "dropped replicate", y = "proportion recovered", fill = NULL) + th

# C: per-tool DEG counts across 2v2 subsets (why the intersection shrinks at n=2)
pC <- summary_df %>% filter(type == "subset_2v2") %>%
  mutate(run = gsub("Huhu.|SUB_", "", run)) %>%
  select(run, DESeq2 = n_deseq2, edgeR = n_edger, `limma-voom` = n_limma, Intersection = n_DEG) %>%
  pivot_longer(-run, names_to = "tool", values_to = "n") %>%
  mutate(tool = factor(tool, levels = c("DESeq2","edgeR","limma-voom","Intersection"))) %>%
  ggplot(aes(reorder(run, n), n, fill = tool)) +
  geom_col(position = position_dodge(.8), width = .75) +
  scale_fill_manual(values = c(DESeq2 = amber, edgeR = "#C9A66B",
                               `limma-voom` = "#B00020", Intersection = ink)) +
  coord_flip() +
  labs(title = "C  DE calls per tool at n=2 (2-vs-2 subsets)",
       subtitle = "limma-voom loses power at n=2, limiting the 3-tool intersection",
       x = NULL, y = "significant transcripts", fill = NULL) +
  th + theme(axis.text.y = element_text(size = 6.5))

# D: read-depth down-sampling
pD <- ds_df %>% group_by(depth) %>%
  summarise(mean = mean(recovery), sd = sd(recovery), .groups = "drop") %>%
  ggplot(aes(depth*100, mean)) +
  geom_ribbon(aes(ymin = pmax(mean-sd,0), ymax = pmin(mean+sd,1)), fill = amber, alpha = .25) +
  geom_line(colour = amber, linewidth = 1) + geom_point(colour = ink, size = 2.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "D  Recovery vs sequencing depth",
       subtitle = "counts thinned to 50-90% (3 reps each)",
       x = "reads retained (%)", y = "DE-set recovery") + th

fig <- (pA | pB) / (pC | pD)
ggsave("sensitivity_out/Figure_Sensitivity.pdf", fig, width = 11, height = 9)
ggsave("sensitivity_out/Figure_Sensitivity.png", fig, width = 11, height = 9, dpi = 200)

# ------------------------------- REPORT --------------------------------------
cat("\n================= SUMMARY =================\n")
cat(sprintf("Reference DE set (3v3 intersection): %d transcripts\n", length(ref_set)))
cat(sprintf("Leave-one-out (n=6):   recovery %.2f-%.2f (median %.2f)\n",
            min(summary_df$recovery[summary_df$type=="LOO"]),
            max(summary_df$recovery[summary_df$type=="LOO"]),
            median(summary_df$recovery[summary_df$type=="LOO"])))
cat(sprintf("2-vs-2 subsets (n=9):  recovery %.2f-%.2f (median %.2f)\n",
            min(summary_df$recovery[summary_df$type=="subset_2v2"]),
            max(summary_df$recovery[summary_df$type=="subset_2v2"]),
            median(summary_df$recovery[summary_df$type=="subset_2v2"])))
cat(sprintf("log2FC Spearman vs full: median %.3f\n", median(summary_df$lfc_spearman, na.rm=TRUE)))
cat(sprintf("Outputs written to: %s/\n", OUTDIR))
cat("==========================================\n")


getwd()                                              # the folder sensitivity_out sits inside
list.files("sensitivity_out", full.names = TRUE)     # the five output files
normalizePath(list.files("sensitivity_out", full.names = TRUE))  # full absolute paths




loo <- summary_df[summary_df$type=="LOO",]
sub <- summary_df[summary_df$type=="subset_2v2",]
cat(sprintf("LOO DE-set recovery %.2f-%.2f (median %.2f)\n", min(loo$recovery), max(loo$recovery), median(loo$recovery)))
cat(sprintf("LOO top-50 recovery %.2f-%.2f\n", min(loo$top_recovery), max(loo$top_recovery)))
cat(sprintf("LOO log2FC Pearson %.2f-%.2f\n", min(loo$lfc_pearson), max(loo$lfc_pearson)))
cat(sprintf("2v2 log2FC Spearman %.2f-%.2f\n", min(sub$lfc_spearman), max(sub$lfc_spearman)))
aggregate(recovery~depth, read.csv("sensitivity_out/downsampling_recovery.csv"), mean)










