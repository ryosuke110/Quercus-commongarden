#!/usr/bin/env Rscript
# Fit a GAM for total dry weight using PCA traits
# Author: Ryosuke Ito
# --- Input files (Dryad) ---
# infile: Imp2507-PC.csv
#
# --- Output files ---
# parametric_outfile: gam_dwtot_pc_parametric.csv
# model_outfile: gam_dwtot_pc_summary.csv

library(mgcv)

### Input ###
infile <- "phenotype-PC.csv"
parametric_outfile <- "gam_dwtot_pc_parametric.csv"
model_outfile <- "gam_dwtot_pc_summary.csv"

### Read data ###
df <- read.csv(infile, check.names = FALSE)

### Prepare data ###
pc_cols <- paste0("PC", 1:10)
need_cols <- c(pc_cols, "HeatShock", "Density", "Dwtot")

df_model <- df[, need_cols]
df_model$HeatShock <- factor(as.integer(df_model$HeatShock))
df_model$Density <- factor(as.integer(df_model$Density))
df_model$Dwtot <- suppressWarnings(as.numeric(df_model$Dwtot))

df_model <- df_model[complete.cases(df_model), ]

### Fit GAM ###
fit <- gam(
  Dwtot ~
    s(PC1) + s(PC2) + s(PC3) + s(PC4) + s(PC5) +
    s(PC6) + s(PC7) + s(PC8) + s(PC9) + s(PC10) +
    HeatShock + Density,
  data = df_model,
  method = "REML"
)

fit_sum <- summary(fit)
gam.check(fit)

### Extract parametric terms ###
parametric_df <- data.frame(
  term = rownames(fit_sum$p.table),
  estimate = fit_sum$p.table[, "Estimate"],
  se = fit_sum$p.table[, "Std. Error"],
  t_value = fit_sum$p.table[, "t value"],
  p_value = fit_sum$p.table[, "Pr(>|t|)"],
  row.names = NULL,
  check.names = FALSE
)

### Extract smooth terms and model summary ###
smooth_df <- data.frame(
  term = rownames(fit_sum$s.table),
  edf = fit_sum$s.table[, "edf"],
  ref_df = fit_sum$s.table[, "Ref.df"],
  F = fit_sum$s.table[, "F"],
  p_value = fit_sum$s.table[, "p-value"],
  n = nrow(df_model),
  r_sq = fit_sum$r.sq,
  dev_expl = fit_sum$dev.expl,
  row.names = NULL,
  check.names = FALSE
)

### Save output ###
write.csv(parametric_df, parametric_outfile, row.names = FALSE)
write.csv(smooth_df, model_outfile, row.names = FALSE)
