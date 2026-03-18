#!/usr/bin/env python3
# Assign each significant SNP to the PC with the largest absolute effect
# Author: Ryosuke Ito

import argparse
import numpy as np
import pandas as pd


def parse_rs(rs: str):
    parts = str(rs).split(":")
    if len(parts) < 3:
        return (np.nan, np.nan, np.nan)
    chr_, pos = parts[1], parts[2]
    return chr_, pos, f"{chr_}:{pos}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input GEMMA result table")
    parser.add_argument("--output", required=True, help="Output trait assignment TSV")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep=",", dtype=str)

    for col in ["q_wald", "beta_1", "beta_2", "beta_3", "beta_4"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    sig = df[df["q_wald"] < args.q_threshold].copy()

    betas = sig[["beta_1", "beta_2", "beta_3", "beta_4"]].abs()
    sig["best_pc"] = "PC" + (betas.values.argmax(axis=1) + 1).astype(str)

    sig[["chr", "pos", "id"]] = pd.DataFrame(
        sig["rs"].apply(parse_rs).tolist(),
        index=sig.index
    )

    sig[["chr", "pos", "id", "best_pc"]].to_csv(
        args.output,
        sep="\t",
        index=False,
        header=False
    )


if __name__ == "__main__":
    main()