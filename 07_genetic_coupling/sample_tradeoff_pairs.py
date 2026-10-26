#!/usr/bin/env python3
# Sample inter-chromosomal SNP pairs for trait-associated and background loci
# Author: Ryosuke Ito

import argparse
import collections
import random


def sample_interchr(pool_a, pool_b, n, max_use=5):
    out = []
    tries = 0
    limit = n * 500
    use_a = collections.Counter()
    use_b = collections.Counter()

    chr_a = [c for c in pool_a if pool_a[c]]
    chr_b = [c for c in pool_b if pool_b[c]]

    if not chr_a or not chr_b:
        return []

    while len(out) < n and tries < limit:
        c1 = random.choice(chr_a)
        c2 = random.choice(chr_b)

        if c1 == c2:
            tries += 1
            continue

        a = random.choice(pool_a[c1])
        b = random.choice(pool_b[c2])

        if a == b or use_a[a] >= max_use or use_b[b] >= max_use:
            tries += 1
            continue

        out.append((a, b))
        use_a[a] += 1
        use_b[b] += 1
        tries += 1

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dp-trait", required=True)
    parser.add_argument("--dp-binned", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--out-cat1", required=True)
    parser.add_argument("--out-cat4", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-use", type=int, default=5)
    parser.add_argument("--target-per-bin", type=int, default=500)
    args = parser.parse_args()

    random.seed(args.seed)

    by_trait_bin_chr = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(list))
    )
    all_by_bin_chr = collections.defaultdict(lambda: collections.defaultdict(list))
    trait_ids = set()

    with open(args.dp_trait) as f:
        for line in f:
            c, pos, sid, af1, af2, dp, bin_, trait = line.rstrip().split("\t")
            by_trait_bin_chr[trait][bin_][c].append(sid)
            trait_ids.add(sid)

    with open(args.dp_binned) as f:
        for line in f:
            c, pos, sid, af1, af2, dp, bin_ = line.rstrip().split("\t")
            if sid not in trait_ids:
                all_by_bin_chr[bin_][c].append(sid)

    bins = []
    x = 0.20
    while x < 0.60 - 1e-9:
        bins.append(f"[{x:.2f},{x+0.05:.2f})")
        x += 0.05
    bins.append("[0.60,0.65]")

    with open(args.pairs) as f:
        pairs = [tuple(line.split()) for line in f if line.strip()]
    pairs = [(a, b) for a, b in pairs if a != b]

    with open(args.out_cat1, "w") as o1, open(args.out_cat4, "w") as o4:
        o1.write("CAT\tPAIR\tBIN\tSNP_A\tSNP_B\n")
        o4.write("CAT\tPAIR\tBIN\tSNP_A\tSNP_B\n")

        for t1, t2 in pairs:
            for bin_ in bins:
                pool_a = by_trait_bin_chr[t1].get(bin_, {})
                pool_b = by_trait_bin_chr[t2].get(bin_, {})

                if not pool_a or not pool_b:
                    continue

                samp1 = sample_interchr(pool_a, pool_b, args.target_per_bin, args.max_use)
                n1 = len(samp1)

                if n1 == 0:
                    continue

                for a, b in samp1:
                    o1.write(f"cat1\t{t1}-{t2}\t{bin_}\t{a}\t{b}\n")

                bg_pool = all_by_bin_chr.get(bin_, {})
                if bg_pool:
                    samp4 = sample_interchr(bg_pool, bg_pool, n1, args.max_use)
                    for a, b in samp4:
                        o4.write(f"cat4\tBG({t1}-{t2})\t{bin_}\t{a}\t{b}\n")


if __name__ == "__main__":
    main()
