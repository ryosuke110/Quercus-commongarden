#!/usr/bin/env Rscript
# Fit a GAM for stomatal conductance with repeated LI-600 measurements
# Author: Ryosuke Ito

library(data.table)
library(tidyr)
library(dplyr)
library(lubridate)
library(mgcv)

### Input ###
infile <- "LI600.txt"

### Read data ###
df <- fread(
  infile,
  sep = "\t",
  na.strings = c("NA", "NaN", "", "na")
)

### Reshape data ###
# Convert paired measurements (_1, _2) into long format
dat_long <- df %>%
  pivot_longer(
    cols = c(
      Time_1, Tleaf_1, VPDleaf_1, PARi_1, gsw_1,
      Time_2, Tleaf_2, VPDleaf_2, PARi_2, gsw_2
    ),
    names_to = c(".value", "rep"),
    names_pattern = "(Time|Tleaf|VPDleaf|PARi|gsw)_(1|2)"
  ) %>%
  mutate(
    Time = suppressWarnings(ymd_hm(Time, tz = "Asia/Tokyo")),
    PlantID = factor(PlantID),
    logPARi = log1p(PARi)
  ) %>%
  filter(
    !is.na(HybridIndex), is.finite(HybridIndex),
    !is.na(gsw), is.finite(gsw),
    !is.na(VPDleaf), is.finite(VPDleaf),
    !is.na(Tleaf), is.finite(Tleaf),
    !is.na(logPARi), is.finite(logPARi)
  )

### Fit GAM ###
# Include major environmental covariates and plant-level random effects
gam_hi_min <- gam(
  gsw ~
    s(HybridIndex, k = 10) +
    s(VPDleaf, k = 10) +
    s(logPARi, k = 10) +
    s(PlantID, bs = "re"),
  data = dat_long,
  method = "REML"
)

summary(gam_hi_min)
AIC(gam_hi_min)
