#!/usr/bin/env Rscript
# Calculate transect distance from geographic coordinates
# Author: Ryosuke Ito

library(data.table)
library(sf)
library(dplyr)

### Input ###
infile   <- "cline.csv"
id_zero  <- "H29301"   # the name of sample at origin of transect distance
crs_geo  <- 4326       # WGS84
crs_prj  <- 6677       # projected CRS; adjust according to study region

### Read data ###
dat <- fread(infile)
stopifnot(all(c("vcfID","HybridIndex","Latitude","Longitude","Elevation") %in% names(dat)))
stopifnot(id_zero %in% dat$vcfID)

### Calculate transect distance ###
pts_sf <- st_as_sf(dat, coords = c("Longitude","Latitude"), crs = crs_geo) |>
  st_transform(crs_prj)
xy <- st_coordinates(pts_sf)

pca <- prcomp(xy, scale. = FALSE)
v1  <- pca$rotation[,1]
origin_xy <- st_coordinates(pts_sf[dat$vcfID == id_zero, ])

signed_dist_m <- as.numeric((xy - matrix(origin_xy, nrow(xy), 2, byrow = TRUE)) %*% v1)

# Flip direction so that distance increases toward the northern side
idx_max_lat <- which.max(dat$Latitude)
vec_to_max  <- xy[idx_max_lat, ] - origin_xy
if (sum(vec_to_max * v1) < 0) signed_dist_m <- -signed_dist_m

dist_m  <- abs(signed_dist_m)
dat$Transect_km <- dist_m / 1000

### Write output ###
outfile <- sub("\\.csv$", "", infile)
outfile <- paste0(outfile, ".with_transect.csv")
fwrite(dat, outfile)
cat("Saved:", outfile, "\n")
