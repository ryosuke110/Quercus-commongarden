#!/usr/bin/env Rscript
# Impute trait data using missForest and evaluate imputation accuracy
# Author: Ryosuke Ito

library(dplyr)
library(missForest)

### Input ###
infile <- "AllData2507.csv"
outfile <- "Imp2507.csv"

trait_cols <- c(
  "HybridIndex", "SampleID", "Height", "Diameter", "CA", "LeafTotal", "DryWeightRoot",
  "DryWeightStem", "DryWeightLeaf", "MeanPredation", "ScanLMA",
  "TensileThickness", "TensileYoungsModulus", "GarminationDate",
  "FilledLeafArea", "LobationIndex", "EFAPC1", "EFAPC2",
  "D.RGR.2208", "D.RGR.2308", "H.RGR.2208", "H.RGR.2308",
  "L.RGR.2208", "L.RGR.2308", "C.RGR.2208", "C.RGR.2309",
  "Density", "HeatShock"
)

### Read data ###
df <- read.csv(infile, check.names = FALSE)

### Prepare data ###
# Keep samples that passed leaf-quality filtering
dat_filt <- df %>%
  filter(Quality == 1) %>%
  select(all_of(trait_cols))

# Exclude identifier and environmental variables from imputation target
dat_imp <- dat_filt %>%
  select(-HybridIndex, -SampleID, -Density, -HeatShock)

### Evaluate imputation accuracy ###
# Use complete cases as pseudo-true data
dat_true <- na.omit(dat_imp)

set.seed(123)
dat_masked <- prodNA(dat_true, noNA = 0.1)

imp_test <- missForest(dat_masked, ntree = 1000, verbose = TRUE)

error_result <- mixError(
  ximp = imp_test$ximp,
  xmis = dat_masked,
  xtrue = dat_true
)

print(error_result)

mask_matrix <- is.na(dat_masked)

imputed_vals <- as.matrix(imp_test$ximp)[mask_matrix]
true_vals <- as.matrix(dat_true)[mask_matrix]

rmse <- sqrt(mean((imputed_vals - true_vals)^2))
nrmse <- rmse / sd(true_vals)
mae <- mean(abs(imputed_vals - true_vals))
correlation <- cor(imputed_vals, true_vals)

cat("NRMSE:", nrmse, "\n")
cat("MAE:", mae, "\n")
cat("Pearson's r:", correlation, "\n")

### Run imputation ###
imp_full <- missForest(dat_imp, ntree = 1000, verbose = TRUE)
dat_completed <- imp_full$ximp

### Save output ###
write.csv(dat_completed, outfile, row.names = FALSE)
