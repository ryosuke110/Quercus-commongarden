#!/usr/bin/env python3
# Reformat fcgene VCF output for downstream analyses
# Author: Ryosuke Ito

import io
import os
import gzip
import argparse
import pandas as pd


def smart_open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def read_vcf(path):
    with smart_open(path, 'rt') as f:
        lines = []
        for line in f:
            if line.startswith('##'):
                continue
            lines.append(line)
        return pd.read_csv(io.StringIO(''.join(lines)), sep='\t')


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input VCF file (.vcf or .vcf.gz)')
parser.add_argument('-o', '--output', required=True, help='Output VCF file (.vcf or .vcf.gz)')
parser.add_argument('-s', '--sample', required=True, help='Sample name file (tab-separated: old_name new_name)')
parser.add_argument('-r', '--refname', required=True, help='Reference name (assembly name)')
parser.add_argument('-n', '--numchr', type=int, default=0, help='Number of contigs')
parser.add_argument('--delsite', action='store_true', help='Delete sites with ALT == "."')
args = parser.parse_args()

### Read VCF ###
# Read VCF header
with smart_open(args.input, 'rt') as f:
    all_lines = f.readlines()

header_lines = [l for l in all_lines if l.startswith('##')]

new_header_lines = []
for line in header_lines:
    if line.startswith('##FILTER=<ID'):
        new_line = line.replace('##FILTER=<ID', '##FORMAT=<ID')
        new_header_lines.append(new_line)
    else:
        new_header_lines.append(line)

for i in range(args.numchr):
    contig_line = f'##contig=<ID=contig{i+1},assembly={args.refname}>\n'
    new_header_lines.append(contig_line)

final_header = ''.join(new_header_lines)

# Read VCF body
dat = read_vcf(args.input)

### Fix VCF ###
# Rename sample columns
sample_df = pd.read_csv(args.sample, sep='\t', header=None, names=['old', 'new'])
for row in sample_df.itertuples(index=False):
    if row.old in dat.columns:
        dat = dat.rename(columns={row.old: row.new})

# Extract CHROM and POS from ID
dat_id = dat['ID'].str.split('_', expand=True)
dat_id.columns = ['CHROM', 'POS']
dat = dat.drop(columns=['#CHROM', 'POS'], errors='ignore')
dat = pd.concat([dat_id, dat], axis=1)
dat = dat.rename(columns={'CHROM': '#CHROM'})

# Optionally remove sites with ALT == "."
if args.delsite:
    dat = dat.query('ALT != "."')

# Replace numeric bases with letters
dat = dat.replace({'REF': {'1': 'C', '2': 'G', '3': 'T', '4': 'A'},
                   'ALT': {'1': 'C', '2': 'G', '3': 'T', '4': 'A'}})

### Write output ###
tmp_header = 'tmp_header.txt'
tmp_data = 'tmp_data.txt'

with open(tmp_header, 'w') as f:
    f.write(final_header)

dat.to_csv(tmp_data, sep='\t', index=False, header=True)

if args.output.endswith('.gz'):
    with open(tmp_header, 'r') as fh, open(tmp_data, 'r') as fd, gzip.open(args.output, 'wt') as fout:
        fout.writelines(fh)
        fout.writelines(fd)
else:
    with open(args.output, 'w') as fout, open(tmp_header, 'r') as fh, open(tmp_data, 'r') as fd:
        fout.writelines(fh)
        fout.writelines(fd)

os.remove(tmp_header)
os.remove(tmp_data)
