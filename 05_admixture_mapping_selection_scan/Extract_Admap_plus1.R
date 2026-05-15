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

### Output ###
outdir <- "adm_plus"
output_files <- list(
  matrix_out = "evidence_matrix.tsv"
  plus1_out = "adm_plus1.tsv"
  plus2_out= "adm_plus2.tsv"
  plus1_pretty_out = "adm_plus1_pretty.tsv"
  plus2_pretty_out = "adm_plus2_pretty.tsv"
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

  dt <- fread(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE
  )

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
    select(gene, n_evidence, evidence)
}

### Build gene sets ###
genes_adm <- read_genes_any(input_files$adm)

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
  genes_adm,
  genes_xp,
  genes_fst,
  genes_dxy,
  genes_pi,
  genes_taj
)))

evidence_df <- data.frame(
  gene = all_genes,
  AdMap = all_genes %in% genes_adm,
  XP = all_genes %in% genes_xp,
  Fst = all_genes %in% genes_fst,
  Dxy = all_genes %in% genes_dxy,
  Pi = all_genes %in% genes_pi,
  TajD = all_genes %in% genes_taj,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

### Count selection evidence ###
evidence_df <- evidence_df %>%
  mutate(
    n_evidence = XP + Fst + Dxy + Pi + TajD
  )

### AdMap + 1 ###
adm_plus1 <- evidence_df %>%
  filter(AdMap, n_evidence >= 1L) %>%
  arrange(desc(n_evidence), gene)

### AdMap + 2 ###
adm_plus2 <- evidence_df %>%
  filter(AdMap, n_evidence >= 2L) %>%
  arrange(desc(n_evidence), gene)

### Save output ###
dir.create(outdir, showWarnings = FALSE)

fwrite(
  evidence_df,
  file = file.path(outdir, output_files$matrix_out),
  sep = "\t"
)

fwrite(
  adm_plus1,
  file = file.path(outdir, output_files$plus1_out),
  sep = "\t"
)

fwrite(
  adm_plus2,
  file = file.path(outdir, output_files$plus2_out),
  sep = "\t"
)

fwrite(
  pretty_flags(adm_plus1),
  file = file.path(outdir, output_files$plus1_pretty_out),
  sep = "\t"
)

fwrite(
  pretty_flags(adm_plus2),
  file = file.path(outdir, output_files$plus2_pretty_out),
  sep = "\t"
)
