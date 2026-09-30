#!/usr/bin/env python3
"""Build gene-level TPM matrix from StringTie -A gene abundance tables.

Usage:
    python3 build_tpm_matrix.py --abund_dir 04.stringtie --output gene_TPM_matrix.tsv
      [--samples sample1,sample2,...]  # restrict to these sample IDs (default: all)
"""
import argparse
import os
import csv
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--abund_dir', required=True,
                    help='Directory containing per-sample *.gene_abund.tab files')
    ap.add_argument('--output', required=True)
    ap.add_argument('--pattern', default='gene_abund.tab',
                    help='File suffix/pattern for abundance tables')
    ap.add_argument('--samples', default=None,
                    help='Comma-separated list of sample IDs to include (default: all found)')
    args = ap.parse_args()

    d = args.abund_dir
    if not os.path.isdir(d):
        sys.exit(f"ERROR: {d} is not a directory")

    include = set(s for s in (args.samples or '').split(',') if s) if args.samples else None

    files = []
    for root, _dirs, names in os.walk(d):
        for fn in names:
            if not fn.endswith(args.pattern):
                continue
            sname = fn.replace('.' + args.pattern, '').replace('_quant', '')
            if include is not None and sname not in include:
                continue
            files.append(os.path.join(root, fn))
    files.sort()
    if not files:
        sys.exit(f"ERROR: no *{args.pattern} files found under {d}")

    print(f"Found {len(files)} abundance tables")

    all_tpm = {}
    for path in files:
        sname = os.path.basename(path).replace('.' + args.pattern, '')
        sname = sname.replace('_quant', '')
        gene_tpm = {}
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gid = row.get('Gene ID', '')
                tpm = float(row.get('TPM', 0) or 0)
                gene_tpm[gid] = tpm
        all_tpm[sname] = gene_tpm
        print(f"  {sname}: {len(gene_tpm)} genes")

    all_genes = sorted(set().union(*[set(t.keys()) for t in all_tpm.values()]))
    print(f"Total unique genes: {len(all_genes)}")

    sample_names = sorted(all_tpm.keys())
    with open(args.output, 'w') as out:
        out.write('gene_id\t' + '\t'.join(sample_names) + '\n')
        for gene in all_genes:
            vals = [str(all_tpm[s].get(gene, 0.0)) for s in sample_names]
            out.write(f"{gene}\t" + '\t'.join(vals) + '\n')

    print(f"TPM matrix written: {args.output} ({len(all_genes)} genes x {len(sample_names)} samples)")


if __name__ == '__main__':
    main()
