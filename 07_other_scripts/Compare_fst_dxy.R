#!/usr/bin/env Rscript
# Compare Fst and Dxy between GWAS outlier windows and background windows
# Author: Ryosuke Ito

library(data.table)

set.seed(42)

### Input ###
outlier_file <- "outliers_gwas_sig.bed"
fst_file <- "fst_10k.csv"
dxy_file <- "dxy_10k.csv"

result_outfile <- "wilcox_outlier_vs_background.tsv"

### Read data ###
outliers <- fread(
  outlier_file,
  col.names = c("CHR", "START", "END")
)

fst <- fread(fst_file)
dxy <- fread(dxy_file)

### Prepare data ###
outliers[, `:=`(
  START = as.integer(START),
  END = as.integer(END)
)]

fst[, `:=`(
  Start = as.integer(Start),
  End = as.integer(End)
)]

dxy[, `:=`(
  Start = as.integer(Start),
  End = as.integer(End)
)]

### Identify outlier windows ###
setkey(outliers, CHR, START, END)
setkey(fst, Chromosome, Start, End)
setkey(dxy, Chromosome, Start, End)

# Convert outlier SNPs to point intervals for overlap with windows
ol_pts <- copy(outliers)[, `:=`(
  snp_start = START,
  snp_end = START
)]
setkey(ol_pts, CHR, snp_start, snp_end)

hit_fst <- foverlaps(
  ol_pts,
  fst,
  by.x = c("CHR", "snp_start", "snp_end"),
  by.y = c("Chromosome", "Start", "End"),
  type = "within",
  nomatch = 0L
)

# Count each window only once even if it contains multiple outlier SNPs
fst_outlier_windows <- unique(
  hit_fst[, .(Chromosome = CHR, Start, End, Fst)],
  by = c("Chromosome", "Start", "End")
)

fst_out_vals <- fst_outlier_windows$Fst
fst_bg_windows <- fst[!fst_outlier_windows, on = .(Chromosome, Start, End)]
fst_bg_vals <- fst_bg_windows$Fst

dxy_out_vals <- dxy[
  fst_outlier_windows,
  on = .(Chromosome, Start, End),
  nomatch = 0L
]$Dxy

dxy_bg_vals <- dxy[
  !fst_outlier_windows,
  on = .(Chromosome, Start, End)
]$Dxy

### Run Wilcoxon tests ###
fst_wilcox <- wilcox.test(
  fst_out_vals,
  fst_bg_vals,
  alternative = "greater"
)

dxy_wilcox <- wilcox.test(
  dxy_out_vals,
  dxy_bg_vals,
  alternative = "greater"
)

### Summarize results ###
result_df <- data.table(
  metric = c("Fst", "Dxy"),
  n_outlier = c(length(fst_out_vals), length(dxy_out_vals)),
  n_background = c(length(fst_bg_vals), length(dxy_bg_vals)),
  median_outlier = c(median(fst_out_vals, na.rm = TRUE), median(dxy_out_vals, na.rm = TRUE)),
  median_background = c(median(fst_bg_vals, na.rm = TRUE), median(dxy_bg_vals, na.rm = TRUE)),
  statistic = c(unname(fst_wilcox$statistic), unname(dxy_wilcox$statistic)),
  p_value = c(fst_wilcox$p.value, dxy_wilcox$p.value)
)

### Save output ###
fwrite(result_df, result_outfile, sep = "\t")
print(result_df)
