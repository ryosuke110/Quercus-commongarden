#!/usr/bin/env python3
# Compare cline models of hybrid index using PyMC
# Author: Ryosuke Ito

import pandas as pd
import numpy as np
import pymc as pm
import arviz as az

### Input ###
admix_file = "cline.with_transect.csv"

### Read data ###
df = pd.read_csv(admix_file)
need = ["HybridIndex", "Elevation", "Transect_km"]
df = df.dropna(subset = need).copy()

elev = df["Elevation"].to_numpy(float)
dist = df["Transect_km"].to_numpy(float)
hi_raw = df["HybridIndex"].to_numpy(float)

### Prepare response ###
# Clip HI values for numerical stability under the Beta likelihood
eps = 1e-6
y = np.clip(hi_raw, eps, 1 - eps)

### Null model ###
with pm.Model() as hi_null_model:
    mu0 = pm.Beta("mu0", alpha = 2, beta = 2)
    kappa0 = pm.Gamma("kappa0", alpha = 2.0, beta = 0.5)

    pm.Beta(
        "HI",
        alpha = mu0 * kappa0,
        beta = (1 - mu0) * kappa0,
        observed = y
    )

    idata_null = pm.sample(
        draws = 20000,
        tune = 20000,
        chains = 4,
        target_accept = 0.995,
        idata_kwargs = {"log_likelihood": True},
        random_seed = 101
    )

### Distance cline model ###
with pm.Model() as hi_dist_model:
    m_d = pm.Normal("center_dist", mu = 1560, sigma = 250)
    w_d = pm.HalfNormal("width_dist", sigma = 5.0)
    s = pm.Normal("slope_sign", 0.0, 1.0)

    eta = s * (dist - m_d) / w_d
    mu = pm.Deterministic("mu", pm.math.sigmoid(eta))
    kappa = pm.Gamma("kappa", alpha = 2.0, beta = 0.5)

    pm.Beta(
        "HI",
        alpha = mu * kappa,
        beta = (1 - mu) * kappa,
        observed = y
    )

    idata_dist = pm.sample(
        draws = 20000,
        tune = 20000,
        chains = 4,
        target_accept = 0.995,
        idata_kwargs = {"log_likelihood": True},
        random_seed = 203
    )

### Elevation cline model ###
with pm.Model() as hi_elev_model:
    m_e = pm.Normal("center_elev", mu = 378.0, sigma = 50.0)
    w_e = pm.HalfNormal("width_elev", sigma = 200.0)

    eta = (elev - m_e) / w_e
    mu = pm.Deterministic("mu", pm.math.sigmoid(eta))
    kappa = pm.Gamma("kappa", alpha = 2.0, beta = 0.5)

    pm.Beta(
        "HI",
        alpha = mu * kappa,
        beta = (1 - mu) * kappa,
        observed = y
    )

    idata_elev = pm.sample(
        draws = 20000,
        tune = 20000,
        chains = 4,
        target_accept = 0.995,
        idata_kwargs = {"log_likelihood": True},
        random_seed = 303
    )

### Compare models ###
cmp = az.compare(
    {"elev": idata_elev, "dist": idata_dist, "null": idata_null},
    ic = "loo",
    method = "BB-pseudo-BMA",
    scale = "deviance"
)

print(cmp)
