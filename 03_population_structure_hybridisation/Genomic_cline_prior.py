#!/usr/bin/env python3
# Estimate a logistic prior for species turnover along temperature
# Author: Ryosuke Ito

import pandas as pd
import numpy as np
import statsmodels.api as sm

### Input ###
infile = "Quercus_climate.csv"

### Read data ###
df = pd.read_csv(infile)
df = df[["species", "temp_annual_mean"]].dropna()
df["species"] = df["species"].astype(int)

# Check that species is binary (0/1) with both classes present
vals = set(df["species"].unique())
if not vals.issubset({0, 1}) or len(vals) < 2:
    raise ValueError(f"species must contain both 0 and 1 (observed: {sorted(vals)})")

### Fit logistic regression ###
X = sm.add_constant(df["temp_annual_mean"].astype(float))
y = df["species"].astype(int)

model = sm.Logit(y, X).fit(disp = False)

alpha = model.params["const"]
beta = model.params["temp_annual_mean"]

cov = model.cov_params().loc[
    ["const", "temp_annual_mean"],
    ["const", "temp_annual_mean"]
].values

### Estimate transition temperature ###
# T* = -alpha / beta
T_star = -alpha / beta

### Calculate 95% CI by the delta method ###
# g(alpha, beta) = -alpha / beta
# dg/dalpha = -1 / beta
# dg/dbeta = alpha / beta^2
grad = np.array([-1.0 / beta, alpha / (beta ** 2)])
var_T = grad @ cov @ grad
se_T = np.sqrt(var_T)

ci_low = T_star - 1.96 * se_T
ci_high = T_star + 1.96 * se_T

### Print results ###
# Logistic regression summary
print(model.summary())

# Estimated transition temperature
print(f"alpha = {alpha:.6f}")
print(f"beta = {beta:.6f}")
print(f"T_star = {T_star:.6f}")
print(f"95% CI = [{ci_low:.6f}, {ci_high:.6f}]")
