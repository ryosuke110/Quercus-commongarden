#!/usr/bin/env Rscript
# Fit GAMs of principal components against hybrid index
# Author: Ryosuke Ito

library(data.table)
library(mgcv)

### Input ###
infile <- "phenotype-PC.csv"
outfile <- "GAM-HI_PCs.csv"
n_pc <- 10

### Read data ###
df <- fread(infile)

### Prepare data ###
df <- df %>%
  transform(
    HybridIndex = suppressWarnings(as.numeric(HybridIndex)),
    Density = as.integer(Density),
    Treatment = factor(as.integer(Treatment))
  )

pc_cols <- grep("^PC\\d+$", names(df), value = TRUE)
pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
pc_cols <- head(pc_cols, n_pc)

need_cols <- c("HybridIndex", "Density", "Treatment", pc_cols)
dat_gam <- df[complete.cases(df[, ..need_cols]), ..need_cols]

### Fit GAMs ###
result_df <- data.frame(
  PC = pc_cols,
  n = NA_integer_,
  HI_edf = NA_real_,
  HI_p = NA_real_,
  HI_p_adj = NA_real_,
  R2 = NA_real_,
  Dev_expl = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(pc_cols)) {
  pc <- pc_cols[i]

  fit <- gam(
    as.formula(paste0(pc, " ~ s(HybridIndex) + Density + Treatment")),
    data = dat_gam,
    method = "REML"
  )

  fit_sum <- summary(fit)
  AIC(fit)

  # Extract the p-value for the smooth term of HybridIndex
  ridx <- grep("^s\\(HybridIndex\\)$", rownames(fit_sum$s.table))

  result_df$n[i] <- nrow(dat_gam)
  result_df$HI_edf[i] <- fit_sum$s.table[ridx, "edf"]
  result_df$HI_p[i] <- fit_sum$s.table[ridx, "p-value"]
  result_df$R2[i] <- fit_sum$r.sq
  result_df$Dev_expl[i] <- fit_sum$dev.expl
}

### Adjust p-values ###
result_df$HI_p_adj <- p.adjust(result_df$HI_p, method = "BH")

### Save output ###
fwrite(result_df, outfile)
