#!/usr/bin/env Rscript
# Fit site-level GAMs for germination and survival
# Author: Ryosuke Ito

library(data.table)
library(dplyr)
library(mgcv)

### Input ###
infile <- "survibval.csv"

### Read data ###
df <- fread(infile, check.names = FALSE)

### Prepare data ###
# Adjust variable types as needed
df <- df %>%
  mutate(
    HybridIndex = as.numeric(HybridIndex),
    Garmination = as.integer(Garmination),
    Survival = as.integer(Survival)
  )

### Summarize by site ###
summary_site <- df %>%
  group_by(SamplingSite) %>%
  summarise(
    n_total = n(),
    n_germinated = sum(Garmination, na.rm = TRUE),
    n_survived = sum(Survival, na.rm = TRUE),
    germ_rate = mean(Garmination, na.rm = TRUE),
    surv_rate = mean(Survival, na.rm = TRUE),
    mean_hindex = mean(HybridIndex, na.rm = TRUE),
    sd_hindex = sd(HybridIndex, na.rm = TRUE),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    .groups = "drop"
  )

# Remove sites with missing mean hybrid index
summary_site_filt <- summary_site %>%
  filter(!is.na(mean_hindex))

### Fit GAMs ###
# Germination rate as a function of site mean hybrid index
gam_germ <- gam(
  cbind(n_germinated, n_total - n_germinated) ~ s(mean_hindex),
  family = binomial(link = "logit"),
  data = summary_site_filt,
  method = "REML"
)

summary(gam_germ)

# Survival rate as a function of site mean hybrid index
# If survival should be conditional on germination, replace n_total with n_germinated
gam_surv <- gam(
  cbind(n_survived, n_total - n_survived) ~ s(mean_hindex),
  family = binomial(link = "logit"),
  data = summary_site_filt,
  method = "REML"
)

summary(gam_surv)
