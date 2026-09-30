#!/usr/bin/env python3
"""Build transcript/protein ID -> gene ID maps from reference GFF/GTF.

Sequences produced by `gffread -w` / `gffread -y` use the transcript or
protein ID from the annotation, e.g. rna-ABC_0001 / prot-ABC_0001.  The
gene_id is found in the GFF/GTF attribute column (or as the Parent feature).

Usage:
  python3 reference_prot_gene_map.py \
      --gff genome.gff --gtf genome.gtf \
      --out protein_to_gene_id.txt \
      --trans_to_gene transcript_to_gene_id.txt

Output:
  protein_id<TAB>gene_id
  transcript_id<TAB>gene_id
"""
import argparse
import re
import sys


def parse_features(path):
    """Yield (feature_type, attribute_dict) for usable annotation lines.

    GFF attributes are ;-separated key=value pairs.  GTF attributes are
    ;-separated key "value"; pairs.  This parser handles both.
    """
    if not path:
        return
    with open(path) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in {'gene', 'mRNA', 'transcript', 'CDS'}:
                continue
            attrs = parts[8]
            fields = {}
            for match in re.finditer(r'([A-Za-z_][A-Za-z0-9_]*)\s*[= ]\s*("([^"]*)"|[^;]+)', attrs):
                key = match.group(1)
                val = match.group(3) if match.group(3) is not None else match.group(2)
                fields[key] = val.strip().rstrip(';')
            yield ftype, fields


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', help='reference GFF (gffread -y/-w source)')
    ap.add_argument('--gtf', help='reference GTF (attribute source, optional)')
    ap.add_argument('--out', required=True,
                    help='protein_id -> gene_id output')
    ap.add_argument('--trans_to_gene',
                    default=None,
                    help='optional transcript_id -> gene_id output')
    args = ap.parse_args()

    if not args.gff and not args.gtf:
        ap.error('need at least --gff or --gtf')

    gene_by_trans = set()          # (transcript_id, gene_id)
    gene_by_prot = set()           # (protein_id, gene_id)
    gene_ids = set()
    transcript_gene_cache = {}

    def add_pair(pairs, key, gid):
        key = key.strip()
        gid = gid.strip()
        if key and gid:
            pairs.add((key, gid))

    # First pass: collect gene-level IDs (GFF) / transcript->gene (GTF).
    if args.gff:
        for ftype, fields in parse_features(args.gff):
            if ftype == 'gene':
                gid = fields.get('ID') or fields.get('gene_id') or fields.get('locus_tag')
                if gid:
                    gene_ids.add(gid)
            elif ftype in {'mRNA', 'transcript'}:
                tid = fields.get('ID') or fields.get('transcript_id')
                parent = fields.get('Parent') or fields.get('gene_id')
                if tid and parent:
                    # Parent can be a comma-separated list; primary is usually first.
                    pid = parent.split(',')[0].strip()
                    gene_ids.add(pid)
                    transcript_gene_cache[tid] = pid
                    add_pair(gene_by_trans, tid, pid)

    if args.gtf:
        for ftype, fields in parse_features(args.gtf):
            gid = fields.get('gene_id', '')
            tid = fields.get('transcript_id', '')
            if tid:
                transcript_gene_cache.setdefault(tid, gid)
                add_pair(gene_by_trans, tid, gid)

    # Second pass over GFF for CDS -> protein_id -> gene_id, using cached
    # transcript->gene relations found above.
    if args.gff:
        for ftype, fields in parse_features(args.gff):
            if ftype != 'CDS':
                continue
            prot = fields.get('protein_id') or fields.get('ID') or ''
            parent = fields.get('Parent') or fields.get('transcript_id') or ''
            tid = parent.split(',')[0].strip()
            gid = transcript_gene_cache.get(tid) or ''
            if prot:
                add_pair(gene_by_prot, prot, gid)

    # gffread -y names sequences by the transcript ID (rna-X) unless the GFF
    # contains a protein_id.  Add rna-/prot- variants so either header works.
    if args.gtf or transcript_gene_cache:
        for tid, gid in list(transcript_gene_cache.items()):
            # gffread -y/-w can keep the plain transcript ID as the FASTA
            # header.  Write the plain transcript ID into the protein lookup
            # too, so downstream ID conversion always works.
            add_pair(gene_by_prot, tid, gid)
            for suffix in ('rna-', 'prot-'):
                add_pair(gene_by_prot, suffix + tid.split(':')[-1], gid)

    def write_pairs(path, pairs, header=None):
        with open(path, 'w') as out:
            seen = set()
            for key, gid in sorted(pairs):
                if (key, gid) in seen:
                    continue
                seen.add((key, gid))
                out.write(f'{key}\t{gid}\n')

    write_pairs(args.out, gene_by_prot)
    if args.trans_to_gene:
        write_pairs(args.trans_to_gene, gene_by_trans)

    print(f'  protein->gene entries: {len(gene_by_prot)}')
    print(f'  transcript->gene entries: {len(gene_by_trans)}')


if __name__ == '__main__':
    main()
