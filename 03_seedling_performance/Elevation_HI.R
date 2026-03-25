#!/usr/bin/env Rscript
# Fit individual- and site-level GAMs for fitness and hybrid index
# Author: Ryosuke Ito
# --- Input files (Dryad) ---
# infile: fitness2507.csv

library(dplyr)
library(mgcv)

### Input ###
infile <- "fitness.csv"

### Read data ###
df <- read.csv(infile, check.names = FALSE)

### Prepare individual-level data ###
dat_ind <- df %>%
  filter(
    !is.na(HybridIndex),
    !is.na(Elevation),
    !is.na(HeatShock)
  ) %>%
  mutate(
    Density = na_if(Density, "na"),
    Density = as.numeric(Density),
	HeatShock = factor(HeatShock)
  ) %>%
  filter(!is.na(Density))

dat_site <- df %>%
  filter(!is.na(HybridIndex), !is.na(Elevation))

### Fit individual-level GAMs ###
gam_ind_hi <- gam(
  HybridIndex ~ s(Elevation) +
    HeatShock + Density,
  data = dat_ind,
  method = "REML"
)

summary(gam_ind_hi)

### Prepare site-level data ###
# Use all individuals with non-missing hybrid index and elevation
# for site-level summaries, regardless of Density availability
site_hi <- dat_site %>%
  group_by(SamplingSite) %>%
  summarise(
    n = n(),
    mean_hindex = mean(HybridIndex, na.rm = TRUE),
    sd_hindex = sd(HybridIndex, na.rm = TRUE),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    sd_elevation = sd(Elevation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_hindex), !is.na(mean_elevation))

### Fit site-level GAM ###
# Site mean hybrid index as a function of site mean elevation
gam_site_hi <- gam(
  mean_hindex ~ s(mean_elevation),
  data = site_hi,
  weights = n,
  method = "REML"
)

summary(gam_site_hi)
