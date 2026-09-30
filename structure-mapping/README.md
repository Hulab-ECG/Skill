# Structure Mapping

Map candidate protein structures to template host enzyme structures using one
of three parallel methods:

1. Foldseek large-scale screening
2. USalign fine structural alignment
3. DEER dense enzyme retrieval (CPU only by default)

## Install

You can provide your own tool paths and conda environments.

- Foldseek: install from the official Foldseek release or provide the path to
  an existing binary.
- USalign: install from the official USalign source/release or provide the path
  to an existing binary.
- DEER: clone [WangJiuming/deer](https://github.com/WangJiuming/deer) and
  create the CPU environment:

```bash
git clone https://github.com/WangJiuming/deer.git
cd deer
conda env create --name deer --file env/env_cpu.yml
conda activate deer
```

If the DEER checkpoint does not auto-download:

```bash
huggingface-cli download cuhkaih/deer --local-dir ./ckpt
```

## Quick run

Foldseek:

```bash
foldseek easy-search   /path/to/candidate_structures/   /path/to/Foldseek_database/   foldseek_hits.m8   tmp_folder   --exhaustive-search 1   -s 7.5   --seq-id-mode 2   --format-output query,target,pident,nident,qcov,qlen,tcov,tlen,evalue,qtmscore,ttmscore,alntmscore,lddt,prob
```

USalign:

```bash
USalign query_structure.pdb   -dir2 candidate_structures/   -split 0   candidate_list   -outfmt 2 > usalign_tmscores.txt
```

DEER CPU:

```bash
CUDA_VISIBLE_DEVICES="" python do_retrieval.py   --template_pdb_dir /path/to/template_pdb/   --database_pdb_dir /path/to/candidate_pdb/   --output_dir /path/to/output
```

Use GPU only if the user explicitly asks for it and the DEER environment was
created with a GPU YAML.

## Outputs

- Foldseek: `foldseek_hits.m8` / `foldseek_hits.tsv`
- USalign: `usalign_results.jsonl`
- DEER: `similarity.csv` sorted by Euclidean embedding distance

See `SKILL.md` for route-selection and filtering guidance.
