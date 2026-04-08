#!/usr/bin/env Rscript
# Run elliptic Fourier analysis and PCA on leaf outlines
# Author: Ryosuke Ito

library(Momocs)
library(dplyr)

### Input ###
infile <- "All_filter.csv"
coord_dir <- "coord"
outfile_pc <- "PCAll-EFA.csv"

### Read metadata ###
dat_main <- read.csv(infile, sep = ",", check.names = FALSE)

### Prepare metadata ###
dat_main_filt <- dat_main %>%
  filter(quality == 1)

### Helper function ###
load_contours <- function(file_names, coord_dir = ".") {
  coord_files <- file.path(coord_dir, file_names)
  import_txt(coord_files, header = FALSE, sep = ",", col.names = c("x", "y"))
}

prepare_shapes <- function(coord_data, slide_dir = "S") {
  outlines <- coo_close(Out(coord_data))
  outlines_centered <- coo_center(outlines)
  outlines_scaled <- Out(
    sapply(outlines_centered$coo, function(x) x / sqrt(coo_area(x)))
  )
  outlines_aligned <- coo_align(outlines_scaled)
  coo_slidedirection(outlines_aligned, slide_dir)
}

run_efa_pca <- function(shape_data, nb_harmonics = 100) {
  efourier_fit <- efourier(shape_data, nb_harmonics, norm = FALSE)
  PCA(efourier_fit)
}

### Quality check ###
# Inspect outlines visually in subsets if needed
qc_files <- dat_main_filt$file_name[1:50]
qc_coord <- load_contours(qc_files, coord_dir = coord_dir)
qc_shape <- prepare_shapes(qc_coord, slide_dir = "S")
panel(qc_shape, names = TRUE)

### Run EFA and PCA ###
coord_data <- load_contours(dat_main_filt$file_name, coord_dir = coord_dir)
shape_data <- prepare_shapes(coord_data, slide_dir = "S")
pca_fit <- run_efa_pca(shape_data, nb_harmonics = 100)

### Plot PCA ###
plot(pca_fit, col = "#066292", cex = 1)

### Save PCA scores ###
write.csv(
  pca_fit$x[, 1:6],
  outfile_pc,
  col.names = TRUE
)

### Remove PC1 outliers ###
dat_main_filt$pc1 <- pca_fit$x[, 1]

pc1_mean <- mean(dat_main_filt$pc1, na.rm = TRUE)
pc1_sd <- sd(dat_main_filt$pc1, na.rm = TRUE)

dat_no_outlier <- dat_main_filt %>%
  filter(
    pc1 >= (pc1_mean - 3 * pc1_sd),
    pc1 <= (pc1_mean + 3 * pc1_sd)
  )

### Re-run EFA and PCA after outlier removal ###
coord_data_filt <- load_contours(dat_no_outlier$file_name, coord_dir = coord_dir)
shape_data_filt <- prepare_shapes(coord_data_filt, slide_dir = "S")
pca_fit_filt <- run_efa_pca(shape_data_filt, nb_harmonics = 100)
