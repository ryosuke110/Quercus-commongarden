#!/usr/bin/env Rscript
# Integrate admixture mapping and selection-scan candidate genes
# Author: Ryosuke Ito

library(data.table)
library(dplyr)
library(stringr)

### Input ###
input_files <- list(
  adm = "outliers_adm_gff.bed",
  xp_pos = "outliers_xpehh_pos_gff.bed",
  xp_neg = "outliers_xpehh_neg_gff.bed",
  fst = "outliers_fst_high_gff.bed",
  dxy = "outliers_dxy_high_gff.bed",
  pi_qmon = "outliers_pi_Qmon_gff.bed",
  pi_qser = "outliers_pi_Qser_gff.bed",
  taj_qmon = "outliers_taj_Qmon_gff.bed",
  taj_qser = "outliers_taj_Qser_gff.bed"
)

outdir <- "adm_plus_alpha"
max_dist_bp <- 5000L
use_dist <- is.finite(max_dist_bp)

weight_vec <- c(
  XP = 1.0,
  Fst = 1.0,
  Dxy = 0.5,
  Pi = 1.0,
  TajD = 1.0
)

### Helper functions ###
extract_gene <- function(x) {
  nm <- str_match(x, "Name=([^;\\t]+)")[, 2]
  ifelse(
    is.na(nm) | nm == "",
    str_match(x, "ID=gene-([^;\\t]+)")[, 2],
    nm
  )
}

read_genes_any <- function(path) {
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = FALSE)

  if (nrow(dt) == 0) {
    return(character(0))
  }

  key_cols <- which(vapply(
    dt,
    function(col) any(str_detect(col, "Name=|ID="), na.rm = TRUE),
    logical(1)
  ))

  if (length(key_cols) == 0) {
    return(character(0))
  }

  attr_vec <- as.character(dt[[tail(key_cols, 1)]])
  gene_vec <- extract_gene(attr_vec)

  unique(gene_vec[!is.na(gene_vec) & nzchar(gene_vec)])
}

read_genes_closest_with_dist <- function(path) {
  dt <- fread(path, sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = FALSE)

  if (nrow(dt) == 0) {
    return(data.frame(gene = character(0), dist = integer(0)))
  }

  key_cols <- which(vapply(
    dt,
    function(col) any(str_detect(col, "Name=|ID="), na.rm = TRUE),
    logical(1)
  ))

  if (length(key_cols) == 0) {
    stop(sprintf("No attribute column was found in: %s", path))
  }

  attr_vec <- as.character(dt[[tail(key_cols, 1)]])
  gene_vec <- extract_gene(attr_vec)
  dist_vec <- as.integer(suppressWarnings(as.numeric(dt[[ncol(dt)]])))

  data.frame(
    gene = gene_vec,
    dist = dist_vec,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(gene) & nzchar(gene))
}

pretty_flags <- function(df) {
  cols <- c("XP", "Fst", "Dxy", "Pi", "TajD")

  df %>%
    mutate(across(all_of(cols), ~ replace(., is.na(.), FALSE))) %>%
    rowwise() %>%
    mutate(
      evidence = {
        hit <- which(c_across(all_of(cols)))
        paste(cols[hit], collapse = ",")
      }
    ) %>%
    ungroup() %>%
    select(gene, min_dist_bp, evidence)
}

### Build gene sets ###
genes_adm_df <- read_genes_closest_with_dist(input_files$adm)
genes_xp <- sort(unique(c(
  read_genes_any(input_files$xp_pos),
  read_genes_any(input_files$xp_neg)
)))
genes_fst <- read_genes_any(input_files$fst)
genes_dxy <- read_genes_any(input_files$dxy)
genes_pi <- sort(unique(c(
  read_genes_any(input_files$pi_qmon),
  read_genes_any(input_files$pi_qser)
)))
genes_taj <- sort(unique(c(
  read_genes_any(input_files$taj_qmon),
  read_genes_any(input_files$taj_qser)
)))

### Build evidence matrix ###
all_genes <- sort(unique(c(
  genes_adm_df$gene,
  genes_xp,
  genes_fst,
  genes_dxy,
  genes_pi,
  genes_taj
)))

evidence_df <- data.frame(
  gene = all_genes,
  AdMap = all_genes %in% genes_adm_df$gene,
  XP = all_genes %in% genes_xp,
  Fst = all_genes %in% genes_fst,
  Dxy = all_genes %in% genes_dxy,
  Pi = all_genes %in% genes_pi,
  TajD = all_genes %in% genes_taj,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

evidence_df <- evidence_df %>%
  left_join(
    genes_adm_df %>%
      group_by(gene) %>%
      summarise(min_dist_bp = min(dist, na.rm = TRUE), .groups = "drop"),
    by = "gene"
  )

### Apply selection rules ###
sel_a <- evidence_df %>%
  mutate(n_plus = XP + Fst + Dxy + Pi + TajD) %>%
  filter(AdMap, n_plus >= 1L, if (use_dist) min_dist_bp <= max_dist_bp else TRUE) %>%
  arrange(desc(n_plus), min_dist_bp)

sel_b <- evidence_df %>%
  mutate(n_plus = XP + Fst + Dxy + Pi + TajD) %>%
  filter(AdMap, n_plus >= 2L, if (use_dist) min_dist_bp <= max_dist_bp else TRUE) %>%
  arrange(desc(n_plus), min_dist_bp)

sel_c <- evidence_df %>%
  mutate(
    score = XP * weight_vec["XP"] +
      Fst * weight_vec["Fst"] +
      Dxy * weight_vec["Dxy"] +
      Pi * weight_vec["Pi"] +
      TajD * weight_vec["TajD"]
  ) %>%
  filter(AdMap, score >= 2.0, if (use_dist) min_dist_bp <= max_dist_bp else TRUE) %>%
  arrange(desc(score), min_dist_bp)

### Save output ###
dir.create(outdir, showWarnings = FALSE)

fwrite(evidence_df, file = file.path(outdir, "evidence_matrix.tsv"), sep = "\t")
fwrite(sel_a %>% arrange(gene), file = file.path(outdir, "adm_plus1.tsv"), sep = "\t")
fwrite(sel_b %>% arrange(gene), file = file.path(outdir, "adm_plus2.tsv"), sep = "\t")
fwrite(sel_c %>% arrange(gene), file = file.path(outdir, "adm_weighted_score.tsv"), sep = "\t")

fwrite(pretty_flags(sel_a), file.path(outdir, "adm_plus1_pretty.tsv"), sep = "\t")
fwrite(pretty_flags(sel_b), file.path(outdir, "adm_plus2_pretty.tsv"), sep = "\t")
fwrite(pretty_flags(sel_c), file.path(outdir, "adm_weighted_pretty.tsv"), sep = "\t")
