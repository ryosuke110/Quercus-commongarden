#!/usr/bin/env Rscript
# Estimate local genomic clines in SNP batches using bgchm
# Author: Ryosuke Ito

library(vcfR)
library(bgchm)

### Input ###
vcf_file <- "aim.CGp.remove.vcf"
hi_file <- "h_est.txt"

# Parent and hybrid sample columns in the VCF object
p0_cols <- c(1, 132, 137, 148, 151, 158)   # Q. serrata parent
p1_cols <- c(1, 18, 21, 22, 45, 49)        # Q. mongolica parent
hyb_cols <- 1:230                          # hybrids + metadata column

### Analysis parameters ###
sdc <- 0.11
sdv <- 0.49
batch_size <- 2000

### Read data ###
vcf_obj <- read.vcfR(vcf_file, verbose = FALSE)

# Divide data into parental and hybrid groups
vcf_p0 <- vcf_obj[, p0_cols]
vcf_p1 <- vcf_obj[, p1_cols]
vcf_hyb <- vcf_obj[, hyb_cols]

### Extract genotype matrices ###
gt_p0 <- extract.gt(vcf_p0, element = "GT", as.numeric = TRUE)
gt_p1 <- extract.gt(vcf_p1, element = "GT", as.numeric = TRUE)
gt_hyb <- extract.gt(vcf_hyb, element = "GT", as.numeric = TRUE)

gtt_p0 <- t(gt_p0)
gtt_p1 <- t(gt_p1)
gtt_hyb <- t(gt_hyb)

### Read prior hybrid index estimates ###
hi_df <- read.table(hi_file, header = TRUE)
hybrid_index <- hi_df[, 1]

### Define SNP batch ###
args <- commandArgs(trailingOnly = TRUE)
batch_id <- as.numeric(args[1])

lower_bound <- (batch_id - 1) * batch_size + 1
upper_bound <- lower_bound + batch_size - 1

# Adjust upper bound if the last batch is incomplete
upper_bound <- min(upper_bound, ncol(gtt_hyb))

locus_idx <- lower_bound:upper_bound
n_loci <- length(locus_idx)

### Subset genotype matrices ###
sub_hyb <- gtt_hyb[, locus_idx, drop = FALSE]
sub_p0 <- gtt_p0[, locus_idx, drop = FALSE]
sub_p1 <- gtt_p1[, locus_idx, drop = FALSE]

### Fit local genomic clines ###
gradient_mat <- matrix(NA, nrow = n_loci, ncol = 3)
center_mat <- matrix(NA, nrow = n_loci, ncol = 3)

for (i in seq_len(n_loci)) {
  fit <- est_genocl(
    Gx = sub_hyb[, i],
    G0 = sub_p0[, i],
    G1 = sub_p1[, i],
    H = hybrid_index,
    model = "genotype",
    ploidy = "diploid",
    hier = FALSE,
    SDc = sdc,
    SDv = sdv
  )

  gradient_mat[i, ] <- fit$gradient
  center_mat[i, ] <- fit$center
}

### Save output ###
outfile_rda <- paste0("clinesOut", batch_id, ".rda")
save(list = ls(), file = outfile_rda)