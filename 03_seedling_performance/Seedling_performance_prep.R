#!/usr/bin/env Rscript
# Summarize site-level fitness data and prepare for GAM analyses
# Author: Ryosuke Ito

# --- Input file (Dryad) ---
# FITNESS_DATA: fitness2507.csv

library(dplyr)

FITNESS_DATA <- "fitness.csv"

## Read input data
df <- read.csv(FITNESS_DATA)

## Check required columns
required_cols <- c("HybridIndex", "Germination", "Alive", "SamplingSite", "Elevation")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

## Standardize column types
df <- df %>%
  rename(Germination = Germination) %>%
  mutate(
    HybridIndex = as.numeric(HybridIndex),
    Germination = as.integer(Germination),
    Alive = as.integer(Alive)
  )

## Summarize by sampling site
summary_site <- df %>%
  group_by(SamplingSite) %>%
  summarise(
    n_total = n(),
    n_germinated = sum(Germination, na.rm = TRUE),
    n_survived = sum(Alive, na.rm = TRUE),
    germ_rate = mean(Germination, na.rm = TRUE),
    mean_Hindex = mean(HybridIndex, na.rm = TRUE),
    sd_Hindex = sd(HybridIndex, na.rm = TRUE),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    .groups = "drop"
  )

## Exclude sites with missing hybrid index
summary_site_complete <- summary_site %>%
  filter(!is.na(mean_Hindex))
