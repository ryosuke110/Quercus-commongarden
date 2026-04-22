#!/usr/bin/env Rscript
# Fit a GAM for photosynthetic rate using LI-6800 measurements
# Author: Ryosuke Ito

library(tidyverse)
library(lubridate)
library(mgcv)

### Input ###
infile <- "LI6800.txt" # "LI6800.txt"

### Read data ###
df6800 <- readr::read_tsv(
  infile,
  show_col_types = FALSE,
  na = c("NA", "NaN", "", "na")
)

### Prepare data ###
dat6800 <- df6800 %>%
  mutate(
    Time = suppressWarnings(ymd_hms(Time, tz = "Asia/Tokyo")),
    logPPFD = log1p(ppfd_umol_m2_s_1)
  ) %>%
  filter(
    !is.na(HybridIndex), is.finite(HybridIndex),
    !is.na(A_Ca420_umol_m2_s_1), is.finite(A_Ca420_umol_m2_s_1),
    !is.na(leaf_vpd_kPa), is.finite(leaf_vpd_kPa),
    !is.na(leaf_temperature_C), is.finite(leaf_temperature_C)
  )

### Inspect data ###
nrow(dat6800)
summary(dat6800$ppfd_umol_m2_s_1)

### Fit GAM ###
# Include major environmental covariates measured during gas exchange
gam_A <- gam(
  A_Ca420_umol_m2_s_1 ~
    s(HybridIndex, k = 10) +
    s(leaf_vpd_kPa, k = 10) +
    s(leaf_temperature_C, k = 10),
  data = dat6800,
  method = "REML"
)

summary(gam_A)
gam.check(gam_A)
