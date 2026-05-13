#!/usr/bin/env Rscript
# Calculate window-based population genetic statistics using GenoPop
# Author: Ryosuke Ito

library(GenoPop)

### Input & output ###
vcf_file <- "all-impb5.CG.cl.vcf.gz"

qmon_file <- "Qmon.txt"
qser_file <- "Qser.txt"
qhyb_file <- "Qhyb.txt"

win_size <- 10000
step_size <- 1
n_threads <- 4

fst_outfile <- "fst_10k.csv"
dxy_outfile <- "dxy_10k.csv"
pi_qmon_outfile <- "pi_Qmon_10k.csv"
pi_qser_outfile <- "pi_Qser_10k.csv"
taj_qmon_outfile <- "tajD_Qmon_10k.csv"
taj_qser_outfile <- "tajD_Qser_10k.csv"

### Read sample sets ###
qmon_ids <- readLines(qmon_file)
qser_ids <- readLines(qser_file)
qhyb_ids <- readLines(qhyb_file)

# Exclude the other parental population and hybrids
exc_qmon <- unique(c(qser_ids, qhyb_ids))
exc_qser <- unique(c(qmon_ids, qhyb_ids))

### Calculate between-population statistics ###
fst_df <- Fst(
  vcf_file,
  qmon_ids,
  qser_ids,
  weighted = TRUE,
  write_log = FALSE,
  window_size = win_size,
  skip_size = step_size,
  threads = n_threads
)

dxy_df <- Dxy(
  vcf_file,
  qmon_ids,
  qser_ids,
  window_size = win_size,
  skip_size = step_size,
  threads = n_threads
)

### Calculate within-population statistics ###
pi_qmon_df <- Pi(
  vcf_file,
  window_size = win_size,
  skip_size = step_size,
  exclude_ind = exc_qmon,
  threads = n_threads
)

pi_qser_df <- Pi(
  vcf_file,
  window_size = win_size,
  skip_size = step_size,
  exclude_ind = exc_qser,
  threads = n_threads
)

taj_qmon_df <- TajimasD(
  vcf_file,
  window_size = win_size,
  skip_size = step_size,
  exclude_ind = exc_qmon,
  threads = n_threads
)

taj_qser_df <- TajimasD(
  vcf_file,
  window_size = win_size,
  skip_size = step_size,
  exclude_ind = exc_qser,
  threads = n_threads
)

### Save output ###
write.csv(fst_df, fst_outfile, row.names = FALSE)
write.csv(dxy_df, dxy_outfile, row.names = FALSE)
write.csv(pi_qmon_df, pi_qmon_outfile, row.names = FALSE)
write.csv(pi_qser_df, pi_qser_outfile, row.names = FALSE)
write.csv(taj_qmon_df, taj_qmon_outfile, row.names = FALSE)
write.csv(taj_qser_df, taj_qser_outfile, row.names = FALSE)
