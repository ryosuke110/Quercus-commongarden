# Quercus-commongarden
## Analytical workflow for the oak hybrid zone study
This repository contains analysis scripts used for the genomic and phenotypic analyses described in:

Ito et al. (in prep.)  
*Molecular evolutionary evidence for coexistence within oak hybrid zones.*

Raw sequencing data are available from DDBJ, and processed datasets are deposited in Dryad.
Quercus_climate.csv is combines occurrence data from the National Survey on the Natural Environment
and climate variables from the NARO Agro-Meteorological Grid Square Data.

---

## Repository structure

The repository is organized by analysis modules corresponding to the main components of the study.

### 01_variant_calling/  
  Scripts for variant calling, filtering, and genotype imputation.

- **Variantcalling_imputation.sh**:
  Pipeline for SNP calling and genotype imputation.
- **fix-fcgene2.py**:
  Utility script used to correct fcGENE output format.
  
### 02_population_structure_hybridisation/  
  Scripts for genomic cline analyses and hybrid index inference.

- **Genomic_cline_prep.R**:
  Preparation of genotype and metadata for cline analysis.
- **Genomic_cline_prior.py**:
  Estimation of priors for genomic cline models.
- **Genomic_cline_pymc.py**:
  Bayesian genomic cline estimation using PyMC.
- **Genomic_cline_hzar_elevation.R**:
  Genomic cline modelling along elevation using the HZAR framework.

### 03_seedling_performance/  
  Analyses of seedling performance traits as a function of hybrid index and environmental variables.

- **Seedling_performance_prep.R**:
  Data preparation for seedling performance analyses.
- **Aarea_HI.R**:
  Leaf area vs hybrid index analysis.
- **DryWeight_HI.R**:
  Seedling biomass vs hybrid index analysis.
- **Elevation_HI.R**:
  Elevation effects on hybrid index.
- **Germination_survival_HI.R**:
  Germination and survival vs hybrid index analyses.
- **Gsw_HI.R**:
  Stomatal conductance vs hybrid index.

### 04_tradeoff/  
Trait correlation and trade-off analyses.

- **DryWeight_traitPCA.R**:
  PCA of biomass allocation traits.
- **HI_traitPCA.R**:
  Relationship between hybrid index and trait PCs.
- **Trait_corplot.R**:
  Correlation analysis among functional traits.
- **Trait_SEM_bestmodel.R**:
  Structural equation modelling of trait relationships.
- **TraitPCA_loading.R**:
  Extraction and visualization of PCA loadings.

### 05_admixture_mapping_selection_scan/
  Genome scans for loci associated with trait variation and signatures of selection.

- **Run_admixture_mapping.sh**:
  Admixture mapping using GEMMA.
- **Calculation_selection_signals.R**:
  Calculation of population genomic statistics (Fst, Dxy, π, Tajima's D).
- **Calculation_XPEHH.R**:
  XP-EHH selection scan.
- **Outlier_detection.R**:
  Detection of genomic outlier loci.
- **Extract_Admap_plus1.R**:
  Integration of admixture mapping hits with other selection signals.

### 06_genetic_scan/  
  Analyses of genetic coupling among trait-associated loci.

- **Run_genetic_coupling.sh**:
  Main pipeline for genetic coupling analysis.
- **assign_trait_maxpc.py**:
  Assigns SNPs to the trait PC with the largest effect size.
- **sample_tradeoff_pairs.py**:
  Sampling of SNP pairs associated with different traits.
- **extract_ld_pairs.py**:
  Extraction of LD values for sampled SNP pairs.

### 07_other_scripts/  
  Additional analyses and utility scripts.

- **Compare_fst_dxy.R**:
  Comparison of Fst and Dxy between GWAS outlier regions and background.
- **Compare_local_cline.R**:
  Comparison of local cline parameters between focal and background loci.
- **Local_cline_fit.R**:
  Estimation of local genomic clines.
- **Run_local_clines_parallel.sh**:
  Parallel execution of local cline estimation.
- **Trait_imputation.R**:
  Imputation of missing trait values.
- **Run_leaf_EFA.R**:
  Exploratory factor analysis of leaf traits.

### README.md  
  Documentation describing the structure and workflow of this repository.

---

## Software and dependencies

The analyses assume that the following software are available in the user environment unless otherwise specified.

* fastp  
* bwa-mem2  
* samtools  
* bcftools  
* bedtools  
* vcftools  
* GATK 3.8
* Picard  
* bamUtil  
* ANGSD  
* Beagle 3.3 & 5.5
* PLINK
* ADMIXTURE  
* GEMMA  
* triangulaR  
* hzar  
* rehh  
* R  
* Python  
* PyMC  
* missForest  
* corrplot  
* piecewiseSEM  
* bgc-hm  
* ShinyGO
