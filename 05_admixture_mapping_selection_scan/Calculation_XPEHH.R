#!/usr/bin/env Rscript
# Calculate XP-EHH between Q. mongolica and Q. serrata across contigs
# Author: Ryosuke Ito

library(vcfR)
library(rehh)

### Input ###
vcf_qmon <- "Qmon.vcf.gz"
vcf_qser <- "Qser.vcf.gz"
outfile <- "xpehh.csv"

chrs <- c(
  "contig1", "contig2", "contig3", "contig4", "contig5", "contig6",
  "contig7", "contig8", "contig9", "contig10", "contig11", "contig12"
)

### Calculate XP-EHH by contig ###
xpe_all <- NULL

for (chr in chrs) {
  message(sprintf("Processing %s ...", chr))

  tryCatch({
    hh_mon <- data2haplohh(
      hap_file = vcf_qmon,
      vcf_reader = "vcfR",
      chr.name = chr,
      polarize_vcf = FALSE,
      min_maf = 0.05,
      min_perc_geno.hap = 100,
      min_perc_geno.mrk = 100,
      remove_multiple_markers = TRUE
    )

    hh_ser <- data2haplohh(
      hap_file = vcf_qser,
      vcf_reader = "vcfR",
      chr.name = chr,
      polarize_vcf = FALSE,
      min_maf = 0.05,
      min_perc_geno.hap = 100,
      min_perc_geno.mrk = 100,
      remove_multiple_markers = TRUE
    )

    sc_mon <- scan_hh(hh_mon, phased = TRUE)
    sc_ser <- scan_hh(hh_ser, phased = TRUE)

    xpe_chr <- ies2xpehh(
      scan_pop1 = sc_mon,
      scan_pop2 = sc_ser,
      popname1 = "Qmon",
      popname2 = "Qser"
    )

    if (!"CHR" %in% names(xpe_chr)) {
      xpe_chr$CHR <- chr
    }

    if (is.null(xpe_all)) {
      xpe_all <- xpe_chr
    } else {
      xpe_all <- rbind(xpe_all, xpe_chr)
    }

  }, error = function(e) {
    message(sprintf("  -> skipped %s (%s)", chr, e$message))
  })
}

### Save output ###
write.csv(xpe_all, outfile, row.names = FALSE)
