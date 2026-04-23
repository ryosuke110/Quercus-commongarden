#!/usr/bin/env Rscript
# Fit a GAM for dry weight against hybrid index
# Author: Ryosuke Ito

library(data.table)
library(dplyr)
library(mgcv)

### Input ###
infile <- "survival.csv"

### Read data ###
df <- fread(infile, check.names = FALSE)

### Prepare data ###
df_clean <- df %>%
  select(HybridIndex, DryWeightTotal) %>%
  mutate(
    HybridIndex = as.numeric(HybridIndex),
    DryWeightTotal = as.numeric(DryWeightTotal)
  ) %>%
  filter(!is.na(HybridIndex), !is.na(DryWeightTotal))

### Fit GAM ###
gam_dw <- gam(
  DryWeightTotal ~ s(HybridIndex),
  data = df_clean,
  method = "REML"
)

summary(gam_dw)
