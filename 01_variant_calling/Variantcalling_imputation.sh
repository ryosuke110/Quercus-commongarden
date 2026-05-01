#!/usr/bin/env bash
# SNP calling and imputation pipeline
# Author: Ryosuke Ito
# The genome and protein sequence data ("Quercus_mongolica_genome.fasta", 
# "Quercus_mongolica_protein.pep", and "Qlobata.v3.0.PCG.prot.fasta") 
# were obtained from publicly available repositories (https://doi.org/10.6084/m9.figshare.11888118.v2; 
# https://valleyoak.ucla.edu/genomic-resources).

set -euo pipefail

# Adjust file paths according to your environment
RAW_REF="Quercus_mongolica_genome.fasata"
REF="Qmon.fa"
REFID="Qmon"
OVERLAPCLIPPED_BAMLIST="all_overlapclipped_bam.txt"
REALIGNED_BAMLIST="all_realigned_bam.txt"
THREADS=12

# Keep only the 12 anchored superscaffolds and rename them to contig01-contig12
seqkit grep \
  -r -p '^Superscaffold[0-9]+$' \
  "${RAW_REF}" \
  | awk '
    /^>/ {
      sub(/^>Superscaffold/, ">contig")
      if ($0 ~ /^>contig[1-9]([^0-9]|$)/) {
        sub(/^>contig/, ">contig0")
      }
    }
    { print }
  ' > "${REF}"

### 1. Build blast database of Q. lobata ###
makeblastdb \
  -in "Qlobata.v3.0.PCG.prot.fasta" \
  -dbtype prot \
  -parse_seqids \
  -out "Qlob_pe" \
  -hash_index

blastp \
  -query "Quercus_mongolica_protein.pep" \
  -db "Qlob_pe" \
  -out "Qmon.p.10.outfmt7" \
  -evalue 1e-10 \
  -outfmt "7 qseqid sseqid length mismatch gapopen qstart qend sstart send evalue" \
  -max_target_seqs 1

### 2. Build reference index ###
samtools faidx "${REF}"

java -jar "picard.jar" CreateSequenceDictionary \
  R="${REF}" \
  O="${REFID}.dict"

bwa-mem2 index "${REF}"

### 3. Read trimming ###
for ID in $(cat "sample_ids.txt"); do

  fastp \
    --cut_right \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    -i "${ID}_R1.fastq.gz" \
    -I "${ID}_R2.fastq.gz" \
    -o "${ID}.R1.fq.gz" \
    -O "${ID}.R2.fq.gz"

done

### 4. Mapping ###
for ID in $(cat "sample_ids.txt"); do

  bwa-mem2 mem \
    -t "${THREADS}" \
    "${REF}" \
    "${ID}.R1.fq.gz" "${ID}.R2.fq.gz" \
    > "${ID}.sam"

  samtools view -bS "${ID}.sam" > "${ID}.bam"

done

### 5. BAM filtering ###
for ID in $(cat "sample_ids.txt"); do

  samtools view -h -q 20 "${ID}.bam" \
    | samtools sort \
    -o "${ID}.sorted.bam"

done

### 6. Add readgroup and remove duplicate ###
for ID in $(cat "sample_ids.txt"); do

  java -jar "picard.jar" AddOrReplaceReadGroups \
    -I "${ID}.sorted.bam" \
    -O "${ID}.readgroup.bam" \
    -RGID "FLOWCELLID" \
    -RGLB "${ID}_library1" \
    -RGPL "DNB" \
    -RGPU "F23FTSAPJT0028" \
    -RGSM "${ID}"

  ## Remove duplicates and print dupstat file
  java -jar "picard.jar" MarkDuplicates \
    -I "${ID}.readgroup.bam" \
    -O "${ID}.dedup.bam" \
    -M "${ID}.dupstat.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --REMOVE_DUPLICATES true

  ## Clip overlapping paired-end reads using bamUtil
  bam clipOverlap \
    --in "${ID}.dedup.bam" \
    --out "${ID}.overlapclipped.bam" \
    --stats

done

### 7. Indel realignment (GATK 3.8) ###
## Create list of potential in-dels
java -jar "GenomeAnalysisTK.jar" \
  -T RealignerTargetCreator \
  -R "${REF}" \
  -I "${OVERLAPCLIPPED_BAMLIST}" \
  -o "indel_realigner.intervals" \
  -drf BadMate

## Run the indel realigner tool
java -jar "GenomeAnalysisTK.jar" \
  -T IndelRealigner \
  -R "${REF}" \
  -I "${OVERLAPCLIPPED_BAMLIST}" \
  -targetIntervals "indel_realigner.intervals" \
  --consensusDeterminationModel USE_READS \
  --nWayOut "_realigned.bam"

### 8. SNP calling using genotype likelihood ###
# Run each chromosome
# -setMinDepth 2X/ind, -setMaxDepth 20X/ind
# -minInd 0.8X OTU
for CHR in {01..12}; do

  angsd \
    -b "${REALIGNED_BAMLIST}" \
    -ref "${REF}" \
    -out "contig${CHR}" \
    -uniqueOnly 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -trim 0 \
    -C 50 \
    -minMapQ 30 \
    -minQ 20 \
    -minInd 808 \
    -setMinDepth 2014 \
    -setMaxDepth 20140 \
    -nThreads 4 \
    -skipTriallelic 1 \
    -doCounts 1 \
    -minMaf 0.05 \
    -SNP_pval 1e-6 \
    -GL 1 \
    -doGlf 2 \
    -doMajorMinor 4 \
    -doMaf 1 \
    -r "contig${CHR}:"

done

### 9. Imputation ###
## Pre-imputation using Beagle3
for CHR in {01..12}; do

  java -jar "beagle-3.3.2/beagle.jar" \
    like="contig${CHR}.beagle.gz" \
    out="contig${CHR}-imp-b3"

  # Convert beagle to vcf
  fcgene \
    --bgl-gprobs "contig${CHR}-imp-b3.contig${CHR}.beagle.gz.gprobs.gz" \
    --oformat vcf \
    --out "contig${CHR}-impb3"

  python "fix-fcgene2.py" \
    -i "contig${CHR}-impb3_vcf.gz" \
    -o "contig${CHR}-impb3-fixed.vcf.gz" \
    -s "all_samples.txt" \
    -r "${REF}" \
    -n 12 \
    --delsite

  ## Filter GP<0.90 sites
  bcftools +setGT \
    "contig${CHR}-impb3-fixed.vcf.gz" \
    -o "contig${CHR}-impb3-filtered.vcf" \
    -- \
    -t q \
    -n . \
    -e 'FORMAT/GP>=0.90'

  ## Post-imputation using beagle5.5
  java -jar "Beagle5.5/beagle.27Feb25.75f.jar" \
    gt="contig${CHR}-impb3-filtered.vcf" \
    out="contig${CHR}-impb5" \
    gp=true

  ## Fix vcf file
  python "fix-fcgene2.py" \
    -i "contig${CHR}-impb5.vcf.gz" \
    -o "contig${CHR}-impb5-fixed.vcf.gz" \
    -s "all_samples2.txt" \
    -r "${REF}" \
    -n 12

done

## Concatenate imputed vcfs
bcftools concat \
  -o "all-impb5.vcf" \
  "contig01-impb5-fixed.vcf.gz" \
  "contig02-impb5-fixed.vcf.gz" \
  "contig03-impb5-fixed.vcf.gz" \
  "contig04-impb5-fixed.vcf.gz" \
  "contig05-impb5-fixed.vcf.gz" \
  "contig06-impb5-fixed.vcf.gz" \
  "contig07-impb5-fixed.vcf.gz" \
  "contig08-impb5-fixed.vcf.gz" \
  "contig09-impb5-fixed.vcf.gz" \
  "contig10-impb5-fixed.vcf.gz" \
  "contig11-impb5-fixed.vcf.gz" \
  "contig12-impb5-fixed.vcf.gz"

### 10. Filter coding region and +-10kb ###
## Prepare Genes BED file
sed -e 's/Superscaffold/contig/g' "Quercus_mongolica_gff3.gff3" > "Qmon.gff3"

grep -P '\tgene\t' "Qmon.gff3" > "Qmon.genes.gff3"

awk 'BEGIN{OFS="\t"} {print $1,$4-1,$5}' "Qmon.genes.gff3" > "Qmon.genes.bed"

bedtools slop \
  -i "Qmon.genes.bed" \
  -g "../../Qmon.rename.fa.fai" \
  -b 10000 \
  > "Qmon.genes.10kb.bed"

bedtools intersect \
  -v \
  -a "all-impb5.vcf" \
  -b "Qmon.genes.10kb.bed" \
  > "tmp.vcf"

sed -ne '/^#/p' "all-impb5.vcf" > "header.vcf"

cat "header.vcf" "tmp.vcf" > "all-impb5.nogenes.vcf"

### 11. LD pruning ###
## Edit chromosome name from contig to number
bcftools annotate \
  --rename-chrs "chr.txt" \
  "all-impb5.nogenes.vcf" \
  > "all-impb5.chr.vcf"

## Preparation for plink file
# Convert VCF to PLINK format
vcftools --vcf "all-impb5.chr.vcf" --plink --out "all-impb5"

## Filtering
plink --file "all-impb5" --indep 50 3 2 --out "all-impb5_ld" --chr-set 12 no-xy no-mt

plink --file "all-impb5" --extract "all-impb5_ld.prune.in" --recode --out "all-impb5_ld_ldpruned"

plink --file "all-impb5_ld_ldpruned" --recode vcf --out "all-impb5_ldpruned"
