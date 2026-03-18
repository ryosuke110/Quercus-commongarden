#!/usr/bin/env python3
# Extract LD values for sampled SNP pairs
# Author: Ryosuke Ito

import argparse
import glob
import pickle


def build_want(id2pos_file, pairs_tsv, out_pkl):
    id2pos = {}
    with open(id2pos_file) as f:
        for line in f:
            chr_, pos, idx = line.rstrip().split("\t")
            id2pos[idx] = (chr_, pos)

    want = set()
    meta = {}

    with open(pairs_tsv) as f:
        next(f)
        for line in f:
            cat, pair, bin_, a, b = line.rstrip().split("\t")
            if a in id2pos and b in id2pos:
                c1, p1 = id2pos[a]
                c2, p2 = id2pos[b]
                want.add((c1, p1, c2, p2))
                want.add((c2, p2, c1, p1))
                meta[(c1, p1, c2, p2)] = (pair, bin_, a, b)
                meta[(c2, p2, c1, p1)] = (pair, bin_, a, b)

    with open(out_pkl, "wb") as f:
        pickle.dump((want, meta), f)


def extract_ld(want_pkl, in_pattern, outfile):
    with open(want_pkl, "rb") as f:
        want, meta = pickle.load(f)

    files = glob.glob(in_pattern)
    if len(files) == 0:
        raise FileNotFoundError(f"No LD files matched: {in_pattern}")

    with open(outfile, "w") as o:
        o.write("PAIR\tBIN\tSNP_A\tSNP_B\tR2\n")

        for inp in files:
            with open(inp) as f:
                header = f.readline().strip().split()
                idx = {name: i for i, name in enumerate(header)}

                i1 = idx.get("CHR_1", 0)
                i2 = idx.get("POS_1", 1)
                i3 = idx.get("CHR_2", 2)
                i4 = idx.get("POS_2", 3)
                i_r = idx.get("R^2", idx.get("R2", 5))

                for line in f:
                    c = line.rstrip().split()
                    key = (c[i1], c[i2], c[i3], c[i4])
                    if key in want:
                        pair, bin_, a, b = meta[key]
                        r2 = c[i_r]
                        o.write(f"{pair}\t{bin_}\t{a}\t{b}\t{r2}\n")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_build = subparsers.add_parser("build")
    p_build.add_argument("--id2pos", required=True)
    p_build.add_argument("--pairs", required=True)
    p_build.add_argument("--output", required=True)

    p_extract = subparsers.add_parser("extract")
    p_extract.add_argument("--want", required=True)
    p_extract.add_argument("--input-pattern", required=True)
    p_extract.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.command == "build":
        build_want(args.id2pos, args.pairs, args.output)
    elif args.command == "extract":
        extract_ld(args.want, args.input_pattern, args.output)


if __name__ == "__main__":
    main()