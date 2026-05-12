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
  select(HybridIndex, DryWeightTotal, Density, Treatment) %>%
  mutate(
    HybridIndex = as.numeric(HybridIndex),
    DryWeightTotal = as.numeric(DryWeightTotal),
    Density = as.numeric(Density),
    Treatment = factor(Treatment)
  ) %>%
  filter(!is.na(HybridIndex), !is.na(DryWeightTotal), !is.na(Density), !is.na(Treatment))

### Fit GAM ###
gam_dw <- gam(
  DryWeightTotal ~ s(HybridIndex) + Density + Treatment,
  data = df_clean,
  method = "REML"
)

summary(gam_dw)
AIC(gam_dw)
