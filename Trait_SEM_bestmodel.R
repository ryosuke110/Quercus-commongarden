#!/usr/bin/env Rscript
# Fit a piecewise SEM for trait trade-offs and biomass allocation
# Author: Ryosuke Ito

suppressPackageStartupMessages({
  library(piecewiseSEM)
  library(dplyr)
  library(readr)
  library(tidyr)
})

### Input ###
infile <- "Imp2507-HI.csv"
coef_outfile <- "coefs_3rd.csv"
r2_outfile <- "r2_3rd.csv"

### Read data ###
df <- readr::read_csv(infile, show_col_types = FALSE)

### Prepare data ###
# Exclude metadata, RGR traits, and selected variables not used in this SEM
dat_sem <- df %>%
  select(-any_of(c("SampleID", "HeatShock", "Density"))) %>%
  select(where(is.numeric)) %>%
  select(-matches("^RGR")) %>%
  select(-any_of(c("Elevation", "Thk", "BBD", "Pred", "PC1ls", "PC2ls")))

### Define variable blocks ###
leaf_fun <- c("LMA", "ILA", "LI", "YM", "CA", "TLN")
arch <- c("Dia", "Ht")
alloc <- c("Dwroot", "Dwstem", "Dwleaf")

### Define structural equations ###
model_formulas <- list(
  # Functional traits
  m_ILA = ILA ~ HybridIndex,
  m_LI  = LI  ~ HybridIndex,
  m_YM  = YM  ~ HybridIndex,
  m_TLN = TLN ~ HybridIndex,

  # Architecture
  m_Dia = Dia ~ CA + LMA + YM + TLN + HybridIndex,
  m_Ht  = Ht  ~ CA + TLN + ILA + HybridIndex,

  # Biomass allocation
  m_Dwroot = Dwroot ~ LMA + CA + TLN + Dia + Ht,
  m_Dwstem = Dwstem ~ ILA + LMA + CA + Dia + Ht + HybridIndex,
  m_Dwleaf = Dwleaf ~ LMA + ILA + CA + TLN + Dia + Ht
)

### Prepare complete-case dataset ###
vars_needed <- unique(unlist(lapply(model_formulas, all.vars)))
vars_needed <- intersect(vars_needed, names(dat_sem))

dat_model <- dat_sem %>%
  select(all_of(vars_needed)) %>%
  mutate(across(everything(), ~ ifelse(is.finite(.), ., NA_real_))) %>%
  drop_na()

cat("N(original) =", nrow(dat_sem), " -> N(complete-case) =", nrow(dat_model), "\n")

### Build linear models ###
mods <- lapply(
  model_formulas,
  function(f) lm(formula = f, data = dat_model, na.action = "na.fail")
)

### Define residual covariances ###
dedup_cov <- function(covs) {
  if (length(covs) == 0) {
    return(covs)
  }

  keys <- vapply(
    covs,
    function(z) {
      v <- all.vars(z)
      paste(sort(v), collapse = "~~")
    },
    ""
  )

  covs[!duplicated(keys)]
}

as_cov <- function(a, b) {
  if (all(c(a, b) %in% names(dat_model))) {
    list(substitute(x %~~% y, list(x = as.name(a), y = as.name(b))))
  } else {
    list()
  }
}

cov_terms <- c(
  as_cov("ILA", "LI"),
  as_cov("CA", "TLN"),
  as_cov("LI", "CA"),
  as_cov("ILA", "CA"),
  as_cov("ILA", "TLN"),
  as_cov("Dwroot", "Dwstem"),
  as_cov("Dwroot", "Dwleaf")
)

cov_terms <- dedup_cov(cov_terms)

### Fit piecewise SEM ###
sem_fit <- do.call(psem, c(mods, cov_terms))

### Extract results ###
coef_df <- coefs(sem_fit)
r2_df <- rsquared(sem_fit)

### Save output ###
readr::write_csv(as.data.frame(coef_df), coef_outfile)
readr::write_csv(as.data.frame(r2_df), r2_outfile)