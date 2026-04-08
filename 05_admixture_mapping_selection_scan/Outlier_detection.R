#!/usr/bin/env Rscript
# Extract outlier loci and windows from admixture mapping and selection scans
# Author: Ryosuke Ito

library(data.table)
library(dplyr)

### Input ###
xpehh_file <- "xpehh_qmon_vs_qser.csv"
admix_file <- "gemma_admixture_mapping.csv"
fst_file <- "fst_10k.csv"
dxy_file <- "dxy_10k.csv"
pi_qmon_file  <- "pi_Qmon_10k.csv"
pi_qser_file  <- "pi_Qser_10k.csv"
taj_qmon_file <- "tajD_Qmon_10k.csv"
taj_qser_file <- "tajD_Qser_10k.csv"

alpha_adm <- 0.05
q_upper <- 0.99
q_lower <- 0.01

### Read data ###
xpehh_df <- fread(xpehh_file)
admix_df <- fread(admix_file)
fst_df <- fread(fst_file)
dxy_df <- fread(dxy_file)
pi_qmon_df  <- fread(pi_qmon_file)
pi_qser_df  <- fread(pi_qser_file)
taj_qmon_df <- fread(taj_qmon_file)
taj_qser_df <- fread(taj_qser_file)

### Helper functions ###
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
label_contig <- function(n) paste0("contig", n)

### Extract XP-EHH outliers ###
thr_xp_pos <- quantile(xpehh_df$xpehh, q_upper, na.rm = TRUE)
thr_xp_neg <- quantile(xpehh_df$xpehh, q_lower, na.rm = TRUE)

xp_pos <- xpehh_df %>%
  filter(is.finite(xpehh), xpehh >= thr_xp_pos) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(pos) - 1L,
    end = as.integer(pos)
  )

xp_neg <- xpehh_df %>%
  filter(is.finite(xpehh), xpehh <= thr_xp_neg) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(pos) - 1L,
    end = as.integer(pos)
  )

### Extract admixture mapping outliers ###
adm_sig <- admix_df %>%
  filter(q_wald < alpha_adm) %>%
  transmute(
    chr = label_contig(to_num(CHR_NUM)),
    start = as.integer(BP) - 1L,
    end = as.integer(BP)
  )

### Extract Fst and dXY outliers ###
thr_fst <- quantile(fst_df$fst, q_upper, na.rm = TRUE)
thr_dxy <- quantile(dxy_df$dxy, q_upper, na.rm = TRUE)

fst_hi <- fst_df %>%
  filter(is.finite(fst), fst >= thr_fst) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

dxy_hi <- dxy_df %>%
  filter(is.finite(dxy), dxy >= thr_dxy) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

### Extract pi outliers ###
thr_pi_qmon <- quantile(pi_qmon_df$pi, q_lower, na.rm = TRUE)
thr_pi_qser <- quantile(pi_qser_df$pi, q_lower, na.rm = TRUE)

pi_qmon <- pi_df %>%
  filter(pop == "Qmon", is.finite(pi), pi <= thr_pi_qmon) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

pi_qser <- pi_df %>%
  filter(pop == "Qser", is.finite(pi), pi <= thr_pi_qser) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

### Extract Tajima's D outliers ###
thr_taj_qmon <- quantile(taj_qmon_df$tajD, q_lower, na.rm = TRUE)
thr_taj_qser <- quantile(taj_qser_df$tajD, q_lower, na.rm = TRUE)

taj_qmon <- taj_df %>%
  filter(pop == "Qmon", is.finite(tajD), tajD <= thr_taj_qmon) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

taj_qser <- taj_df %>%
  filter(pop == "Qser", is.finite(tajD), tajD <= thr_taj_qser) %>%
  transmute(
    chr = label_contig(CHR_NUM),
    start = as.integer(start) - 1L,
    end = as.integer(end)
  )

### Save output ###
fwrite(xp_pos,   "outliers_xpehh_pos.bed",  sep = "\t", col.names = FALSE)
fwrite(xp_neg,   "outliers_xpehh_neg.bed",  sep = "\t", col.names = FALSE)
fwrite(adm_sig,  "outliers_adm_sig.bed",    sep = "\t", col.names = FALSE)
fwrite(fst_hi,   "outliers_fst_high.bed",   sep = "\t", col.names = FALSE)
fwrite(dxy_hi,   "outliers_dxy_high.bed",   sep = "\t", col.names = FALSE)
fwrite(pi_qmon,  "outliers_pi_Qmon.bed",    sep = "\t", col.names = FALSE)
fwrite(pi_qser,  "outliers_pi_Qser.bed",    sep = "\t", col.names = FALSE)
fwrite(taj_qmon, "outliers_taj_Qmon.bed",   sep = "\t", col.names = FALSE)
fwrite(taj_qser, "outliers_taj_Qser.bed",   sep = "\t", col.names = FALSE)
