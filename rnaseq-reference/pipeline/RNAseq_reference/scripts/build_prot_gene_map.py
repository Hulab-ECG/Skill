#!/usr/bin/env python3
import re
import sys

merged_gtf = sys.argv[1]
out_file = sys.argv[2]
mapping = {}
with open(merged_gtf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 9:
            continue
        if fields[2] != 'transcript':
            continue
        attrs = fields[8]
        g = re.search(r'gene_id "([^"]+)"', attrs)
        t = re.search(r'transcript_id "([^"]+)"', attrs)
        if g and t:
            mapping[t.group(1)] = g.group(1)
with open(out_file, 'w') as out:
    for tid, gid in sorted(mapping.items()):
        out.write(f'{tid}\t{gid}\n')
print(f'  transcript->gene mapping: {len(mapping)} entries')
