#!/usr/bin/env Rscript
# Run PCA on imputed trait data and export PC scores
# Author: Ryosuke Ito
#
# --- Input file (Dryad) ---
# infile: Imp2507-HI.csv
#
# --- Output file ---
# score_outfile: Imp2507-PC.csv
# loading_outfile: Imp2507-PC-loadings.csv

library(dplyr)
library(tibble)

### Input ###
infile <- "phenotype-HI.csv"
score_outfile <- "phenotype-PC.csv"
loading_outfile <- "phenoype-PC-loadings.csv"

trait_cols <- c(
  "Ht", "Dia", "CA", "TLN", "Dwroot", "Dwstem", "Dwleaf", "RSratio",
  "Pred", "LMA", "Thk", "YM", "BBD", "ILA", "LI", "PC1ls", "PC2ls",
  "RGRd2208", "RGRd2308", "RGRht2208", "RGRht2308",
  "RGRtln2208", "RGRtln2308", "RGRca2208", "RGRca2309"
)

meta_cols <- c("SampleID", "HybridIndex", "Dwtot", "HeatShock", "Density")

### Read data ###
df <- read.csv(infile, check.names = FALSE)

### Prepare PCA input ###
dat_pca <- df %>%
  select(all_of(trait_cols))

### Run PCA ###
pca_res <- prcomp(dat_pca, center = TRUE, scale. = TRUE)

### Extract loadings ###
loading_df <- as_tibble(pca_res$rotation, rownames = "Trait")

### Extract PC scores ###
pc_score_df <- as_tibble(pca_res$x[, 1:10])
colnames(pc_score_df) <- paste0("PC", 1:10)

### Combine metadata and PC scores ###
df_pca_summary <- df %>%
  select(all_of(meta_cols)) %>%
  bind_cols(pc_score_df)

### Save output ###
write.csv(df_pca_summary, score_outfile, row.names = FALSE)
write.csv(loading_df, loading_outfile, row.names = FALSE)
