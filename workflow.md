# Workflow

This document describes the analytical workflow associated with:

Ito et al. (in prep.)
*Molecular evolutionary evidence for coexistence within oak hybrid zones.*

Processed datasets are available through Dryad.

Raw sequencing reads have been deposited in DDBJ under BioProject PRJDB40283.

---

# Demo datasets

The processed datasets deposited in this repository can be used as demonstration datasets for reproducing representative analyses described in the associated GitHub repository.

The repository contains processed datasets corresponding to all major analytical modules, including:

* population structure and genomic cline analyses
* seedling performance analyses
* trait trade-off analyses
* admixture mapping and selection scans
* genetic coupling analyses

Detailed instructions for running representative analyses, expected outputs, and typical runtimes are provided in the corresponding sections below.

These processed datasets allow reproduction of representative analyses without reprocessing raw sequencing data from DDBJ.

---

# 1. Variant calling and genotype imputation

This module generates the filtered and imputed SNP dataset used in all downstream analyses.

## Step 1.1 Variant calling and genotype imputation

Script:

01_variant_calling/Variantcalling_imputation.sh

Input:

* Raw sequencing reads (DDBJ)

Output:

* all-CG.vcf.gz

Description:

Reads are quality filtered, aligned to the reference genome, subjected to variant calling, and subsequently imputed to generate the final SNP dataset.

---

# 2. Population structure and genomic cline analyses

This module estimates hybrid indices and genomic cline parameters.

## Step 2.1 Data preparation

Script:

02_population_structure_hybridisation/Genomic_cline_prep.R

Input:

* cline.csv
* sample metadata

Output:

* cline.with_transect.csv

Description:

Converts sample coordinates into one-dimensional transect positions and prepares datasets for genomic cline analyses.

---

## Step 2.2 Prior estimation

Script:

02_population_structure_hybridisation/Genomic_cline_prior.py

Input:

* cline.with_transect.csv

Output:

* prior parameter estimates

Description:

Estimation of prior distributions for Bayesian genomic cline analyses.

---

## Step 2.3 Bayesian genomic cline analysis

Script:

02_population_structure_hybridisation/Genomic_cline_pymc.py

Input:

* cline.with_transect.csv
* prior parameter estimates

Output:

* posterior distributions of cline parameters

Description:

Bayesian estimation of genomic cline parameters using PyMC.

---

## Step 2.4 Elevational genomic cline analysis

Script:

02_population_structure_hybridisation/Genomic_cline_hzar_elevation.R

Input:

* cline.with_transect.csv

Output:

* posterior distributions of cline parameters

Description:

Genomic cline modelling along elevational gradients using the HZAR framework.

---

# 3. Seedling performance analyses

This module evaluates relationships among hybrid index, seedling performance, and physiological traits.

## Step 3.1 Data preparation

Script:

03_seedling_performance/Seedling_performance_prep.R

Input:

* survival.csv

Output:

* survival.by_site.csv

Description:

Preparation of phenotypic datasets used in downstream analyses.

---

## Step 3.2 Germination and survival analyses

Script:

03_seedling_performance/Germination_survival_HI.R

Input:

* survival.csv

Output:

* model summaries
* statistical results

Description:

Analysis of germination and survival as a function of hybrid index.

---

## Step 3.3 Biomass analyses

Script:

03_seedling_performance/DryWeight_HI.R

Input:

* survival.csv

Output:

* model summaries
* statistical results

Description:

Analysis of seedling biomass variation along the hybrid index gradient.

---

## Step 3.4 Photosynthetic rate analyses

Script:

03_seedling_performance/Aarea_HI.R

Input:

* LI6800.txt

Output:

* model summaries
* statistical results

Description:

Analysis of photosynthetic performance measured using the LI-6800 system.

---

## Step 3.5 Stomatal conductance analyses

Script:

03_seedling_performance/Gsw_HI.R

Input:

* LI600.txt

Output:

* model summaries
* statistical results

Description:

Analysis of stomatal conductance variation along the hybrid index gradient.

---

# 4. Trait trade-off analyses

This module evaluates correlations among functional traits, trait principal components, and structural relationships among traits.

## Step 4.1 Trait correlation analysis

Script:

04_tradeoff/Trait_corplot.R

Input:

* phenotype-HI.csv

Output:

* corplot.pdf

Description:

Correlation analysis among functional traits and generation of a trait correlation matrix.

---

## Step 4.2 PCA loading extraction and visualization

Script:

04_tradeoff/TraitPCA_loading.R

Input:

* phenotype-HI.csv

Output:

* phenotype-PC.csv
* phenotype-PC-loadings.csv

Description:

Principal component analysis of trait data and extraction of loading matrices.

---

## Step 4.3 Biomass trait analyses

Script:

04_tradeoff/DryWeight_traitPCA.R

Input:

* phenotype-PC.csv

Output:

* gam_dwtot_pc_parametric.csv
* gam_dwtot_pc_summary.csv

Description:

Analysis of relationships between biomass allocation and trait principal components.

---

## Step 4.4 Trait PCA and hybrid index associations

Script:

04_tradeoff/HI_traitPCA.R

Input:

* phenotype-PC.csv

Output:

* GAM-HI_PCs.csv

Description:

Analysis of relationships between hybrid index and trait principal components.

---

## Step 4.5 Structural equation modelling

Script:

04_tradeoff/Trait_SEM_bestmodel.R

Input:

* phenotype-HI.csv

Output:

* coefs.csv
* r2.csv

Description:

Structural equation modelling of relationships among hybrid index, functional traits, and seedling performance.

---

# 5. Admixture mapping and selection scans

This module identifies loci associated with phenotypic variation and evaluates signatures of divergent selection.

## Step 5.1 Admixture mapping

Script:

05_admixture_mapping_selection_scan/Run_admixture_mapping.sh

Input:

* all-CG.vcf.gz
* all_phenotype-PC.csv
* randomeff_bimbam.txt
* quercus_snpannotation.txt

Output:

* gemma_admixture_mapping.csv

Description:

Admixture mapping using GEMMA to identify SNPs associated with phenotypic trait variation.

---

## Step 5.2 Calculation of population genomic statistics

Script:

05_admixture_mapping_selection_scan/Calculation_selection_signals.R

Input:

* all-CG.vcf.gz
* Qmon.txt
* Qser.txt
* Qhyb.txt

Output:

* fst_10k.csv
* dxy_10k.csv
* pi_Qmon_10k.csv
* pi_Qser_10k.csv
* tajD_Qmon_10k.csv
* tajD_Qser_10k.csv

Description:

Calculation of window-based population genomic statistics, including Fst, Dxy, nucleotide diversity, and Tajima's D.

---

## Step 5.3 XP-EHH analysis

Script:

05_admixture_mapping_selection_scan/Calculation_XPEHH.R

Input:

* Qmon.vcf.gz
* Qser.vcf.gz

Output:

* xpehh.csv

Description:

XP-EHH analysis to detect genomic regions showing evidence of divergent selection between parental lineages.

---

## Step 5.4 Outlier detection

Script:

05_admixture_mapping_selection_scan/Outlier_detection.R

Input:

* gemma_admixture_mapping.csv
* xpehh.csv
* fst_10k.csv
* dxy_10k.csv
* pi_Qmon_10k.csv
* pi_Qser_10k.csv
* tajD_Qmon_10k.csv
* tajD_Qser_10k.csv

Output:

* outlier BED files

Description:

Detection of genomic outlier regions associated with admixture mapping and population genomic selection statistics.

---

## Step 5.5 Integration of admixture mapping and selection signals

Script:

05_admixture_mapping_selection_scan/Extract_Admap_plus1.R

Input:

* outlier BED files

Output:

* adm_plus/evidence_matrix.tsv
* adm_plus/adm_plus1.tsv
* adm_plus/adm_plus2.tsv
* adm_plus/adm_plus1_pretty.tsv
* adm_plus/adm_plus2_pretty.tsv

Description:

Integration of admixture mapping hits with selection-scan outliers and SNP annotation information.

---

# 6. Genetic coupling analyses

This module evaluates genetic coupling among loci associated with different trait axes.

## Main pipeline

Script:

06_genetic_coupling/Run_genetic_coupling.sh

Input:

* admixture mapping results
* genotype data
* trait-PC associations

Output:

* ld_cat1_HZ.tsv
* ld_cat4_HZ.tsv
* ld_cat1_Qmon.tsv
* ld_cat4_Qmon.tsv
* ld_cat1_Qser.tsv
* ld_cat4_Qser.tsv
* pairs_cat1_tradeoff.tsv
* pairs_cat4_bg.tsv
* pair_counts.tsv

Description:

Evaluates genetic coupling among trait-associated loci by comparing inter-chromosomal linkage disequilibrium (ICLD) among trade-off-associated loci against background expectations.

### Internal scripts

assign_trait_maxpc.py

: Assigns each SNP to the trait principal component with the largest effect size.

sample_tradeoff_pairs.py

: Generates SNP pairs associated with different trait axes and background SNP pairs.

extract_ld_pairs.py

: Extracts inter-chromosomal LD values for sampled SNP pairs.

---

# 7. Additional analyses and utility scripts

This directory contains supplementary analyses and utility scripts used in the manuscript. The scripts are largely independent and are not intended to be executed as a single analytical pipeline.

## Local genomic cline fitting

Script:

07_other_scripts/Local_cline_fit.R

Input:

* aim.CGp.remove.vcf
* h_est.txt

Output:

* clinesOut.rda

Description:

Estimation of local genomic cline parameters for ancestry-informative markers.

---

## Parallel execution of local cline analyses

Script:

07_other_scripts/Run_local_clines_parallel.sh

Input:

* aim.CGp.remove.vcf

Output:

* clinesOut.rda

Description:

Parallel execution of local genomic cline analyses.

---

## Comparison of local cline parameters

Script:

07_other_scripts/Compare_local_cline.R

Input:

* mvlmm-query_with_bgchm-0b.csv
* bgchm_merged.csv

Output:

* fg_bg_angle_difference_summary.tsv

Description:

Compares local cline parameters between focal loci and background loci.

---

## Comparison of Fst and Dxy

Script:

07_other_scripts/Compare_fst_dxy.R

Input:

* outliers_gwas_sig.bed
* fst_10k.csv
* dxy_10k.csv

Output:

* wilcox_outlier_vs_background.tsv

Description:

Compares population differentiation and sequence divergence between GWAS outlier regions and genomic background windows.

---

## Trait imputation

Script:

07_other_scripts/Trait_imputation.R

Input:

* all_phenotype.csv

Output:

* imputed phenotype dataset

Description:

Imputation of missing trait values prior to downstream trait analyses.

---

## Leaf shape analysis

Script:

07_other_scripts/Run_leaf_EFA.R

Input:

* coord.zip
* All_filter.csv

Output:

* leaf shape principal components
* elliptic Fourier analysis outputs

Description:

Elliptic Fourier analysis of leaf outline coordinates after filtering low-quality leaf images.

---

# Notes

The scripts are organized according to the main analytical modules of the manuscript. Some scripts require processed input files deposited in Dryad, whereas raw read processing starts from sequencing data deposited in DDBJ. Intermediate temporary files generated during script execution are not exhaustively listed here. These processed datasets allow reproduction of representative analyses on a standard desktop computer without reprocessing raw sequencing data from DDBJ.
