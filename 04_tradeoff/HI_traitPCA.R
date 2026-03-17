#!/usr/bin/env Rscript
# Fit GAMs of principal components against hybrid index
# Author: Ryosuke Ito

library(mgcv)

### Input ###
infile <- "Imp2507-PC.csv"
outfile <- "GAM-HI_PCs_mgcv.csv"
n_pc <- 10

### Read data ###
df <- read.csv(infile, check.names = FALSE)

### Prepare data ###
df <- df %>%
  transform(
    HybridIndex = suppressWarnings(as.numeric(HybridIndex)),
    Density = factor(as.integer(Density)),
    HeatShock = factor(as.integer(HeatShock))
  )

pc_cols <- grep("^PC\\d+$", names(df), value = TRUE)
pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
pc_cols <- head(pc_cols, n_pc)

need_cols <- c("HybridIndex", "Density", "HeatShock", pc_cols)
dat_gam <- df[complete.cases(df[, need_cols]), need_cols]

### Fit GAMs ###
result_df <- data.frame(
  PC = pc_cols,
  HI_p = NA_real_,
  HI_p_adj = NA_real_,
  R2 = NA_real_,
  Dev_expl = NA_real_,
  EDF_total = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(pc_cols)) {
  pc <- pc_cols[i]

  fit <- gam(
    as.formula(paste0(pc, " ~ s(HybridIndex) + Density + HeatShock")),
    data = dat_gam,
    method = "REML"
  )

  fit_sum <- summary(fit)

  # Extract the p-value for the smooth term of HybridIndex
  ridx <- grep("^s\\(HybridIndex\\)$", rownames(fit_sum$s.table))

  result_df$HI_p[i] <- fit_sum$s.table[ridx, "p-value"]
  result_df$R2[i] <- fit_sum$r.sq
  result_df$Dev_expl[i] <- fit_sum$dev.expl
  result_df$EDF_total[i] <- sum(fit_sum$edf)
}

### Adjust p-values ###
result_df$HI_p_adj <- p.adjust(result_df$HI_p, method = "BH")

### Save output ###
write.csv(result_df, outfile, row.names = FALSE)
