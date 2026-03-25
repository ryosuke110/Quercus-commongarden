#!/usr/bin/env Rscript
# Infer clines along elevation and distance gradients using hzar
# Author: Ryosuke Ito
# --- Input files (Dryad) ---
# infile: sample.cline.with_transect.csv

library(hzar)

### Input ###
infile <- "cline.with_transect.csv"

### Read data ###
dat <- read.csv(infile, check.names = FALSE)

# Define site IDs from geographic coordinates
# Adjust rounding precision if needed
dat$site_id <- with(
  dat,
  paste0(sprintf("%.5f", Latitude), "_", sprintf("%.5f", Longitude))
)

### Prepare elevation-based observations ###
# Mean elevation for each site
site_ele <- tapply(dat$Elevation, dat$site_id, mean, na.rm = TRUE)

site_ele_vec <- hzar.mapSiteDist(
  siteID = names(site_ele),
  distance = as.numeric(site_ele)
)

obs_ele <- hzar.doNormalData1DRaw(
  site.dist = site_ele_vec,
  traitSite = dat$site_id,
  traitValue = dat$HybridIndex
)

### Build elevation cline model ###
x <- as.numeric(unname(site_ele_vec))
rng <- range(x, finite = TRUE)
pad <- 0.20 * diff(rng)

model_ele <- hzar.makeCline1DNormal(obs_ele)
model_ele <- hzar.model.addBoxReq(model_ele, rng[1] - pad, rng[2] + pad)

### Set initial parameters for elevation model ###
v_by_site <- tapply(dat$HybridIndex, dat$site_id, var, na.rm = TRUE)
var0 <- median(v_by_site, na.rm = TRUE)
var0 <- ifelse(!is.finite(var0) || var0 <= 1e-4, 1e-3, var0)

muL0 <- max(0, quantile(dat$HybridIndex, 0.05, na.rm = TRUE))
muR0 <- min(1, quantile(dat$HybridIndex, 0.95, na.rm = TRUE))
center0 <- median(x)
width0 <- diff(rng) / 3

pt <- model_ele$parameterTypes

pt$muL$val <- muL0
pt$muL$w <- 0.03
pt$muR$val <- muR0
pt$muR$w <- 0.03
pt$center$val <- center0
pt$center$w <- 0.003 * diff(rng)
pt$width$val <- width0
pt$width$w <- 0.001 * diff(rng)

pt$varL$val <- var0
pt$varL$w <- 0.2 * var0
pt$varR$val <- var0
pt$varR$w <- 0.2 * var0
pt$varH$val <- var0
pt$varH$w <- 0.2 * var0

pt$muL$lower <- 0
pt$muL$upper <- 1
pt$muR$lower <- 0
pt$muR$upper <- 1
pt$varL$lower <- 1e-6
pt$varR$lower <- 1e-6
pt$varH$lower <- 1e-8

model_ele$parameterTypes <- pt

### Fit elevation cline model ###
fitR_ele <- hzar.first.fitRequest.gC(model_ele, obs_ele, verbose = FALSE)
fitR_ele$mcmcParam$chainLength <- 1e7
fitR_ele$mcmcParam$burnin <- 1e6
fitR_ele$mcmcParam$thin <- 1000

fit_ele <- hzar.doFit(fitR_ele)

fitR2_ele <- hzar.next.fitRequest(fit_ele)
fitR2_ele <- hzar.multiFitRequest(fitR2_ele, each = 4, baseSeed = NULL)
runs_ele <- hzar.doChain.multi(fitR2_ele, doPar = FALSE, inOrder = FALSE, count = 4)

### Summarize elevation model ###
summary(
  do.call(mcmc.list, lapply(runs_ele, function(x) hzar.mcmc.bindLL(x[[3]])))
)

init_dgs_ele <- list(
  normal = hzar.dataGroup.add(fit_ele)
)

oDG_ele <- hzar.make.obsDataGroup(init_dgs_ele)
oDG_ele <- hzar.copyModelLabels(init_dgs_ele, oDG_ele)
oDG_ele <- hzar.make.obsDataGroup(lapply(runs_ele, hzar.dataGroup.add), oDG_ele)

AICc_ele <- hzar.AICc.hzar.obsDataGroup(oDG_ele)
print(AICc_ele)

best_ele <- oDG_ele$data.groups[[rownames(AICc_ele)[which.min(AICc_ele$AICc)]]]
print(hzar.get.ML.cline(best_ele))

### Prepare distance-based observations ###
# Mean transect distance for each site
site_dist <- tapply(dat$transect_km, dat$site_id, mean, na.rm = TRUE)

site_dist_vec <- hzar.mapSiteDist(
  siteID = names(site_dist),
  distance = as.numeric(site_dist)
)

obs_dist <- hzar.doNormalData1DRaw(
  site.dist = site_dist_vec,
  traitSite = dat$site_id,
  traitValue = dat$HybridIndex
)

### Build distance cline model ###
x <- as.numeric(unname(site_dist_vec))
rng <- range(x, finite = TRUE)
pad <- 0.20 * diff(rng)

model_dist <- hzar.makeCline1DNormal(obs_dist)
model_dist <- hzar.model.addBoxReq(model_dist, rng[1] - pad, rng[2] + pad)

### Set initial parameters for distance model ###
v_by_site <- tapply(dat$HybridIndex, dat$site_id, var, na.rm = TRUE)
var0 <- median(v_by_site, na.rm = TRUE)
var0 <- ifelse(!is.finite(var0) || var0 <= 1e-4, 1e-3, var0)

muL0 <- max(0, quantile(dat$HybridIndex, 0.05, na.rm = TRUE))
muR0 <- min(1, quantile(dat$HybridIndex, 0.95, na.rm = TRUE))
center0 <- median(x)
width0 <- diff(rng) / 3

pt <- model_dist$parameterTypes

pt$muL$val <- muL0
pt$muL$w <- 0.03
pt$muR$val <- muR0
pt$muR$w <- 0.03
pt$center$val <- center0
pt$center$w <- 0.003 * diff(rng)
pt$width$val <- width0
pt$width$w <- 0.001 * diff(rng)

pt$varL$val <- var0
pt$varL$w <- 0.2 * var0
pt$varR$val <- var0
pt$varR$w <- 0.2 * var0
pt$varH$val <- var0
pt$varH$w <- 0.2 * var0

pt$muL$lower <- 0
pt$muL$upper <- 1
pt$muR$lower <- 0
pt$muR$upper <- 1
pt$varL$lower <- 1e-6
pt$varR$lower <- 1e-6
pt$varH$lower <- 1e-8

model_dist$parameterTypes <- pt

### Fit distance cline model ###
fitR_dist <- hzar.first.fitRequest.gC(model_dist, obs_dist, verbose = FALSE)
fitR_dist$mcmcParam$chainLength <- 1e7
fitR_dist$mcmcParam$burnin <- 1e6
fitR_dist$mcmcParam$thin <- 1000

fit_dist <- hzar.doFit(fitR_dist)

fitR2_dist <- hzar.next.fitRequest(fit_dist)
fitR2_dist <- hzar.multiFitRequest(fitR2_dist, each = 4, baseSeed = NULL)
runs_dist <- hzar.doChain.multi(fitR2_dist, doPar = FALSE, inOrder = FALSE, count = 4)

### Summarize distance model ###
summary(
  do.call(mcmc.list, lapply(runs_dist, function(x) hzar.mcmc.bindLL(x[[3]])))
)

init_dgs_dist <- list(
  normal = hzar.dataGroup.add(fit_dist)
)

oDG_dist <- hzar.make.obsDataGroup(init_dgs_dist)
oDG_dist <- hzar.copyModelLabels(init_dgs_dist, oDG_dist)
oDG_dist <- hzar.make.obsDataGroup(lapply(runs_dist, hzar.dataGroup.add), oDG_dist)

AICc_dist <- hzar.AICc.hzar.obsDataGroup(oDG_dist)
print(AICc_dist)

best_dist <- oDG_dist$data.groups[[rownames(AICc_dist)[which.min(AICc_dist$AICc)]]]
print(hzar.get.ML.cline(best_dist))
