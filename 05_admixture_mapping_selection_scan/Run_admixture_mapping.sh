#!/usr/bin/env bash
# Prepare input files and run admixture mapping with GEMMA
# Author: Ryosuke Ito

set -euo pipefail

### Input ###
# The input VCF is assumed to have already excluded low-quality or unwanted individuals
VCF_IN="all-impb5.CG.vcf"
PREFIX="Quercus_CG"
BED_PREFIX="Quercus_QC"
KEEP_FILE="quercus_with_phenos.txt"

BIMBAM_OUT="Quercus_bimbam_dosage.txt"
BIMBAM_OUT_AD="Quercus_bimbam_dosage_ad.txt"
SNPANN_OUT="oak_snpannotation.txt"
PHENO_FILE="Imp2507-PC-formatted.csv"
KIN_OUT="relmatrix_Q"
GEMMA_OUT="oak_lmm_PCs_2508"

### Prepare PLINK files for ADMIXTURE ###
plink \
  --allow-extra-chr \
  --vcf "${VCF_IN}" \
  --make-bed \
  --out "${PREFIX}"

admixture "${PREFIX}.bed" 2

### Prepare genotype files for GEMMA ###
plink \
  --allow-extra-chr \
  --vcf "${VCF_IN}" \
  --dog \
  --geno 0.0 \
  --make-bed \
  --out "${BED_PREFIX}"

plink \
  --bfile "${BED_PREFIX}" \
  --dog \
  --geno 0.0 \
  --keep "${KEEP_FILE}" \
  --recode \
  --out "${BED_PREFIX}"

### Convert genotype data to BIMBAM dosage format ###
qctool \
  -g "${BED_PREFIX}.bed" \
  -ofiletype bimbam_dosage \
  -og "${BIMBAM_OUT}"

### Add covariate information if needed ###
cat "${BIMBAM_OUT}" randomeff_bimbam.txt > "${BIMBAM_OUT_AD}"

### Calculate relatedness matrix ###
gemma \
  -g "${BIMBAM_OUT_AD}" \
  -notsnp \
  -p "${PHENO_FILE}" \
  -a "${SNPANN_OUT}" \
  -gk 1 \
  -o "${KIN_OUT}"

### Run association analysis ###
gemma \
  -g "${BIMBAM_OUT_AD}" \
  -p "${PHENO_FILE}" \
  -a "${SNPANN_OUT}" \
  -k "output/${KIN_OUT}.cXX.txt" \
  -notsnp \
  -lmm 4 \
  -n 1 2 3 4 \
  -o "${GEMMA_OUT}"
