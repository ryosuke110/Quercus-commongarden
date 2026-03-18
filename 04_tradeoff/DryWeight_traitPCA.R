#!/usr/bin/env Rscript
# Fit a GAM for total dry weight using PCA traits
# Author: Ryosuke Ito

library(mgcv)

### Input ###
infile <- "Imp2507-PC.csv"
rds_outfile <- "gam_dwtot_pc.rds"

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

summary(fit)
gam.check(fit)

### Save model ###
saveRDS(fit, rds_outfile)
