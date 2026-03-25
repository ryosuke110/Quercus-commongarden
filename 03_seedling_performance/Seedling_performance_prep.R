#!/usr/bin/env Rscript
# Fit a GAM for photosynthetic rate using LI-6800 measurements
# Author: Ryosuke Ito
# --- Input files (Dryad) ---
# infile: fitness2507.csv

library(dplyr)
library(mgcv)

infile = "fitness.csv"

## データ読み込み
df <- read.csv(infile)

## 型を整える（HybridIndexを数値に, Garminationを0/1に）
df <- df %>%
  mutate(
    HybridIndex = as.numeric(HybridIndex),
    Garmination = as.integer(Garmination)
  )

## サイトごとの集計
summary_site <- df %>%
  group_by(SamplingSite) %>%
  summarise(
    n_total       = n(),
    n_germinated  = sum(Garmination, na.rm = TRUE),
    n_survived  = sum(Alive, na.rm = TRUE),
    germ_rate     = mean(Garmination, na.rm = TRUE),
    mean_Hindex   = mean(HybridIndex, na.rm = TRUE),
    sd_Hindex     = sd(HybridIndex, na.rm = TRUE),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    .groups = "drop"
  )

## HI が欠損のサイトはモデルから除外
summary_site_isna <- summary_site %>%
  filter(!is.na(mean_Hindex))
