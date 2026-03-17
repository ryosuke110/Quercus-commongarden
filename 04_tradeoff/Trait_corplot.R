#!/usr/bin/env Rscript
# Create a Spearman correlation plot for imputed trait data
# Author: Ryosuke Ito

library(tidyverse)
library(Hmisc)
library(corrplot)

### Input ###
infile <- "Imp2507-HI.csv"
plot_outfile <- "corplot_spearman_ellipse_sigdot.pdf"

drop_cols <- c("SampleID", "HybridIndex", "Elevation", "HeatShock", "Density")
use_abs <- FALSE

### Read data ###
dat <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

### Prepare trait matrix ###
# Exclude identifier and environmental columns, then keep numeric traits only
dat_traits <- dat %>%
  select(-any_of(drop_cols)) %>%
  select(where(is.numeric))

### Calculate correlations ###
rc <- Hmisc::rcorr(as.matrix(dat_traits), type = "spearman")
cor_mat <- rc$r
p_mat <- rc$P

### Reorder traits by hierarchical clustering ###
dist_mat <- if (use_abs) {
  as.dist(1 - abs(cor_mat))
} else {
  as.dist(1 - cor_mat)
}

hc <- hclust(dist_mat, method = "complete")
ord <- hc$order

cor_mat_ord <- cor_mat[ord, ord]
p_mat_ord <- p_mat[ord, ord]

### Define helper function ###
# Add significance dots to the upper triangle only
add_sig_dots <- function(p_mat, alpha = 0.05, cex = 1.2, col = "black") {
  n <- nrow(p_mat)
  idx <- which(p_mat <= alpha & row(p_mat) < col(p_mat), arr.ind = TRUE)

  if (nrow(idx) > 0) {
    xs <- idx[, 2]
    ys <- n - idx[, 1] + 1
    text(xs, ys, labels = "・", cex = cex, col = col)
  }
}

### Plot correlation matrix ###
corrplot(
  cor_mat_ord,
  method = "ellipse",
  type = "upper",
  order = "original",
  tl.col = "black",
  tl.cex = 0.7,
  addgrid.col = "grey90",
  sig.level = 0.05,
  insig = "blank"
)

add_sig_dots(p_mat_ord, alpha = 0.05, cex = 1.2, col = "black")

### Save plot ###
pdf(plot_outfile, width = 9, height = 9)

corrplot(
  cor_mat_ord,
  method = "ellipse",
  type = "upper",
  order = "original",
  tl.col = "black",
  tl.cex = 0.7,
  addgrid.col = "grey90",
  sig.level = 0.05,
  insig = "blank"
)

add_sig_dots(p_mat_ord, alpha = 0.05, cex = 1.2, col = "black")

dev.off()

### Optional dendrogram ###
plot(
  as.dendrogram(hc),
  main = if (use_abs) "Dendrogram (1 - |Spearman r|)" else "Dendrogram (1 - Spearman r)",
  ylab = if (use_abs) "Height (1 - |r|)" else "Height (1 - r)"
)
