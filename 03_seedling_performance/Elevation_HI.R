#!/usr/bin/env Rscript
# Fit individual- and site-level GAMs for fitness and hybrid index
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
  mutate(
    HybridIndex = suppressWarnings(as.numeric(as.character(HybridIndex))),
    Elevation = suppressWarnings(as.numeric(as.character(Elevation))),
    Density = suppressWarnings(as.numeric(as.character(Density))),
    Treatment = factor(Treatment),
    SamplingSite = factor(SamplingSite)
  )

### Prepare individual-level data ###
dat_ind <- df_clean %>%
  filter(
    !is.na(HybridIndex),
    !is.na(Elevation),
	!is.na(Density),
    !is.na(Treatment)
  )

dat_site <- df_clean %>%
  filter(
	  !is.na(HybridIndex), 
	  !is.na(Elevation)
  )

### Fit individual-level GAMs ###
gam_ind_hi <- gam(
  HybridIndex ~ s(Elevation) +
    Treatment + Density,
  data = dat_ind,
  method = "REML"
)

summary(gam_ind_hi)
AIC(gam_ind_hi)

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
AIC(gam_site_hi)
