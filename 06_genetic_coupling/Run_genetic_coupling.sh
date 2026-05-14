#!/usr/bin/env bash
# Compare inter-chromosomal LD between trait-associated and background loci
# Author: Ryosuke Ito

set -euo pipefail

### Input ###
VCF_IN="all-CG.vcf.gz"
OUTDIR="mvlmm"

KEEP_HYB="Qhyb.txt"
KEEP_P1="Qmon.txt"
KEEP_P2="Qser.txt"

TRAIT_CSV="mvlmm_pwald_FDR05_hits.tsv"

SEED=42
MAX_USE=5
TARGET_PER_BIN=500

LD_MODE="hap"      # geno | hap
MIN_R2="0"

### Analysis parameters ###
DELTA_P_MIN="0.20"
DELTA_P_MAX="0.65"
DELTA_P_STEP="0.05"

### Python helper scripts ###
PY_ASSIGN="assign_trait_maxpc.py"
PY_SAMPLE="sample_tradeoff_pairs.py"
PY_EXTRACT="extract_ld_pairs.py"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

### (0) Prepare filtered VCF and SNP ID map ###
# The input VCF is assumed to contain the focal set of individuals.
# SNP IDs are reset to CHROM:POS after filtering.

log "(0) prepare filtered VCF"

vcftools \
  --gzvcf "../${VCF_IN}" \
  --maf 0.1 \
  --max-missing 0.9 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --stdout \
  | bgzip > filtered.vcf.gz

tabix -p vcf filtered.vcf.gz

bcftools annotate -x ID filtered.vcf.gz \
  | bcftools annotate --set-id '%CHROM:%POS' -Oz -o filtered.id.vcf.gz

tabix -p vcf filtered.id.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\n' filtered.id.vcf.gz > id2pos.tsv
awk '{print $3"\t"$1"\t"$2}' id2pos.tsv > id2pos.map

### (1) Calculate parental allele frequencies and bin SNPs by delta-p ###
log "(1) calculate parental allele frequencies and delta-p bins"

bcftools view -S "${KEEP_P1}" -Oz -o P1.vcf.gz filtered.id.vcf.gz
bcftools index -f P1.vcf.gz
bcftools +fill-tags P1.vcf.gz -Oz -o P1.af.vcf.gz -- -t AF
bcftools index -f P1.af.vcf.gz

bcftools view -S "${KEEP_P2}" -Oz -o P2.vcf.gz filtered.id.vcf.gz
bcftools index -f P2.vcf.gz
bcftools +fill-tags P2.vcf.gz -Oz -o P2.af.vcf.gz -- -t AF
bcftools index -f P2.af.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' P1.af.vcf.gz > P1.af.tsv
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' P2.af.vcf.gz > P2.af.tsv

join -t $'\t' -1 3 -2 3 \
  <(sort -t $'\t' -k3,3 P1.af.tsv) \
  <(sort -t $'\t' -k3,3 P2.af.tsv) \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      # join output: ID CHR1 POS1 AF1 CHR2 POS2 AF2
      id=$1; chr1=$2; pos1=$3; af1=$4; chr2=$5; pos2=$6; af2=$7;
      dp=af1-af2; if (dp<0) dp=-dp;
      if (dp < 0.20 || dp > 0.65) next;
      eps=1e-9;
      if (dp >= 0.60 - eps) {
          lower=0.60; upper=0.65; right="]";
      } else {
          lower = 0.20 + 0.05 * int((dp-0.20)/0.05 + eps);
          upper = lower + 0.05; right=")";
      }
      bin=sprintf("[%.2f,%.2f%s", lower, upper, right);
      print chr1, pos1, id, af1, af2, dp, bin
  }' > dp_binned.tsv

### (2) Assign each significant SNP to the PC with the largest absolute effect ###
log "(2) assign trait SNPs to PCs"

python3 "${PY_ASSIGN}" \
  --input "../${TRAIT_CSV}" \
  --output "trait_maxpc.tsv" \
  --q-threshold 0.05

### (3) Merge delta-p bins and trait labels ###
log "(3) merge delta-p bins and trait labels"

join -t $'\t' -1 3 -2 3 \
  <(sort -t $'\t' -k3,3 dp_binned.tsv) \
  <(sort -t $'\t' -k3,3 trait_maxpc.tsv) \
  | awk -F'\t' -v OFS='\t' '{
      # join output: ID + dp table + trait table
      id=$1; chr=$2; pos=$3; af1=$4; af2=$5; dp=$6; bin=$7; trait=$10;
      print chr, pos, id, af1, af2, dp, bin, trait
  }' > dp_trait.tidy.tsv

awk -F'\t' -v OFS='\t' '
  {k=$8 FS $7; c[k]++}
  END{
    print "TRAIT\tBIN\tN_SNP";
    for(k in c) print k, c[k];
  }' dp_trait.tidy.tsv \
  | sort -t $'\t' -k1,1 -k2,2 \
  > counts_trait_bin.tsv

### (4) Generate all trait pairs ###
log "(4) generate trait pairs"

python3 - << 'PY'
import itertools

seen = set()
with open("dp_trait.tidy.tsv") as f:
    for line in f:
        t = line.rstrip().split("\t")[-1]
        seen.add(t)

uniq = sorted(seen, key=lambda x: int(x[2:]) if x.startswith("PC") else 99)
pairs = list(itertools.combinations(uniq, 2))

with open("pairs_tradeoff.txt", "w") as o:
    for a, b in pairs:
        o.write(f"{a} {b}\n")
PY

### (5) Sample inter-chromosomal SNP pairs ###
log "(5) sample trait-associated and background SNP pairs"

python3 "${PY_SAMPLE}" \
  --dp-trait "dp_trait.tidy.tsv" \
  --dp-binned "dp_binned.tsv" \
  --pairs "pairs_tradeoff.txt" \
  --out-cat1 "pairs_cat1_tradeoff.tsv" \
  --out-cat4 "pairs_cat4_bg.tsv" \
  --seed "${SEED}" \
  --max-use "${MAX_USE}" \
  --target-per-bin "${TARGET_PER_BIN}"

### (6) Prepare position sets for LD calculation ###
log "(6) prepare SNP position sets"

tail -n +2 pairs_cat1_tradeoff.tsv | awk '{print $4; print $5}' | sort -u > cat1.ids
awk 'NR==FNR{m[$1]=$2"\t"$3; next} {if($1 in m) print m[$1]}' id2pos.map cat1.ids \
  | sort -u > cat1.union.pos
wc -l cat1.union.pos

tail -n +2 pairs_cat4_bg.tsv | awk '{print $4; print $5}' | sort -u > cat4.ids
awk 'NR==FNR{m[$1]=$2"\t"$3; next} {if($1 in m) print m[$1]}' id2pos.map cat4.ids \
  | sort -u > cat4.union.pos
wc -l cat4.union.pos

### (7) Compute inter-chromosomal LD ###
log "(7) compute inter-chromosomal LD"

ld_flag="--interchrom-geno-r2"
if [ "${LD_MODE}" = "hap" ]; then
  ld_flag="--interchrom-hap-r2"
fi

run_ld() {
  local keep=$1
  local pos=$2
  local prefix=$3

  vcftools \
    --gzvcf filtered.id.vcf.gz \
    --keep "${keep}" \
    --positions "${pos}" \
    ${ld_flag} \
    --ld-window 999999999 \
    --ld-window-bp 999999999 \
    --min-r2 "${MIN_R2}" \
    --out "${prefix}"
}

# Hybrid zone
run_ld "${KEEP_HYB}" cat1.union.pos HZ_cat1
run_ld "${KEEP_HYB}" cat4.union.pos HZ_cat4

# Parental populations
run_ld "${KEEP_P1}" cat1.union.pos Qmon_cat1
run_ld "${KEEP_P1}" cat4.union.pos Qmon_cat4
run_ld "${KEEP_P2}" cat1.union.pos Qser_cat1
run_ld "${KEEP_P2}" cat4.union.pos Qser_cat4

### (8) Extract only sampled LD pairs ###
log "(8) extract sampled LD pairs"

python3 "${PY_EXTRACT}" build \
  --id2pos "id2pos.tsv" \
  --pairs "pairs_cat1_tradeoff.tsv" \
  --output "want_cat1.pkl"

python3 "${PY_EXTRACT}" build \
  --id2pos "id2pos.tsv" \
  --pairs "pairs_cat4_bg.tsv" \
  --output "want_cat4.pkl"

# Hybrid zone
python3 "${PY_EXTRACT}" extract \
  --want "want_cat1.pkl" \
  --input-pattern "HZ_cat1.interchrom.*.ld" \
  --output "ld_cat1_HZ.tsv"

python3 "${PY_EXTRACT}" extract \
  --want "want_cat4.pkl" \
  --input-pattern "HZ_cat4.interchrom.*.ld" \
  --output "ld_cat4_HZ.tsv"

# Q. mongolica
python3 "${PY_EXTRACT}" extract \
  --want "want_cat1.pkl" \
  --input-pattern "Qmon_cat1.interchrom.*.ld" \
  --output "ld_cat1_Qmon.tsv"

python3 "${PY_EXTRACT}" extract \
  --want "want_cat4.pkl" \
  --input-pattern "Qmon_cat4.interchrom.*.ld" \
  --output "ld_cat4_Qmon.tsv"

# Q. serrata
python3 "${PY_EXTRACT}" extract \
  --want "want_cat1.pkl" \
  --input-pattern "Qser_cat1.interchrom.*.ld" \
  --output "ld_cat1_Qser.tsv"

python3 "${PY_EXTRACT}" extract \
  --want "want_cat4.pkl" \
  --input-pattern "Qser_cat4.interchrom.*.ld" \
  --output "ld_cat4_Qser.tsv"

### (9) Summarize sampled pair counts ###
log "(9) summarize sampled pair counts"

{
  echo -e "CAT\tPAIR\tBIN\tN_pairs"
  awk 'NR>1{k=$2 FS $3; c[k]++} END{for(k in c) print "cat1\t" k "\t" c[k]}' pairs_cat1_tradeoff.tsv
  awk 'NR>1{k=$2 FS $3; c[k]++} END{for(k in c) print "cat4\t" k "\t" c[k]}' pairs_cat4_bg.tsv
} | sort -t $'\t' -k1,1 -k2,2 -k3,3 > pair_counts.tsv

log "DONE"
