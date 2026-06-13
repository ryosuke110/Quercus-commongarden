# Quercus-commongarden

Analysis scripts associated with:

Ito et al. (in prep.)
*Molecular evolutionary evidence for coexistence within oak hybrid zones.*

This repository contains scripts used for population genomic, phenotypic, and statistical analyses described in the manuscript.

Processed datasets are available through Dryad.

Raw sequencing reads have been deposited in DDBJ under BioProject PRJDB40283 and will be released upon publication.

Detailed analytical procedures are described in workflow.md.

Quercus_climate.csv combines occurrence data from the National Survey on the Natural Environment and climate variables from the NARO Agro-Meteorological Grid Square Data.

---

## Repository structure

- 01_variant_calling/
- 02_population_structure_hybridisation/
- 03_seedling_performance/
- 04_tradeoff/
- 05_admixture_mapping_selection_scan/
- 06_genetic_coupling/
- 07_other_scripts/

A detailed description of each module is provided below.

---

#### 01_variant_calling/  
  Scripts for variant calling, filtering, and genotype imputation.
- **Variantcalling_imputation.sh**:
  Pipeline for SNP calling and genotype imputation.
- **fix-fcgene2.py**:
  Utility script used to correct fcGENE output format.
  
#### 02_population_structure_hybridisation/  
  Scripts for genomic cline analyses and hybrid index inference.
- **Genomic_cline_prep.R**:
  Preparation of genotype and metadata for cline analysis.
- **Genomic_cline_prior.py**:
  Estimation of priors for genomic cline models.
- **Genomic_cline_pymc.py**:
  Bayesian genomic cline estimation using PyMC.
- **Genomic_cline_hzar_elevation.R**:
  Genomic cline modelling along elevation using the HZAR framework.

#### 03_seedling_performance/  
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

#### 04_tradeoff/  
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

#### 05_admixture_mapping_selection_scan/
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

#### 06_genetic_coupling/  
  Analyses of genetic coupling among trait-associated loci.
- **Run_genetic_coupling.sh**:
  Main pipeline for genetic coupling analysis.
- **assign_trait_maxpc.py**:
  Assigns SNPs to the trait PC with the largest effect size.
- **sample_tradeoff_pairs.py**:
  Sampling of SNP pairs associated with different traits.
- **extract_ld_pairs.py**:
  Extraction of LD values for sampled SNP pairs.

#### 07_other_scripts/  
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

#### README.md  
  Documentation describing the structure and workflow of this repository.

#### workflow.md
  Detailed analytical workflow, including representative commands, input datasets, and expected outputs for each analysis.

---

## Installation

Clone this repository:
```
git clone https://github.com/ryosuke110/Quercus-commongarden.git
cd Quercus-commongarden
```
Typical installation time is less than 1 minute on a standard desktop computer.

## Software requirements

The analyses assume that the following software are available in the user environment unless otherwise specified.

* R (v4.4.0)
* Python (v3.12)
* blastp (v.2.16.0)
* fastp  (v1.0.0)
* bwa-mem2  (v2.2.1)
* Picard  (v3.4.0)
* bamUtil (v1.0.15)
* GATK (v3.8.1)
* ANGSD  (v0.94)
* Beagle (v3.3.1 & v5.5)
* fcGENE (v1.0.7)
* samtools  (v1.17)
* bcftools (v1.22) 
* bedtools (v2.18)
* PLINK (v1.9 & v2.0)
* vcftools (v0.1.17)
* ADMIXTURE (v1.3.0)
* triangulaR (v1.14)
* hzar (v0.2)
* PyMC (v5.26)
* mgcv (v1.9.4)
* missForest (v1.6.0)
* Hmisc (v5.2.0)
* corrplot (v0.95)
* piecewiseSEM (v2.3.1)
* GEMMA (v0.98)
* rehh (v3.2.2)
* GenoPop (v0.9)
* ShinyGO (v0.85)
* bgc-hm (version not specified)
* ggplot2 (4.0.0)

---

## Analytical workflow

The workflow consists of six major analytical modules:

1. Variant calling and genotype imputation
2. Population structure and genomic cline analyses
3. Seedling performance analyses
4. Trait trade-off analyses
5. Admixture mapping and selection scans
6. Genetic coupling analyses

A detailed description of the workflow, including representative commands, input datasets, and expected outputs, is provided in workflow.md.

---

## License

This repository is distributed under the MIT License.
