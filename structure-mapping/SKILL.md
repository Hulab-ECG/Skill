---
name: structure-mapping
description: Map candidate protein structures to template host enzyme structures using Foldseek large-scale screening, USalign fine structural alignment, or DEER dense enzyme retrieval. Use when candidate structures must be narrowed to high-confidence hits.
---

# Structure Mapping

Given candidate protein structures and template host enzyme structures, this
skill produces a shortlist of structurally similar candidates. Step 1 (Foldseek),
Step 2 (USalign), and Step 3 (DEER) are parallel options; pick one route based
on scale, available tools, and the user's preference.

## Inputs

- Candidate structure directory containing PDB or CIF files.
- Template host enzyme structure directory or a list of template PDB files.

## Outputs

- Foldseek route: `foldseek_hits.m8`, `foldseek_hits.tsv`, and a Foldseek-based
  candidate list filtered by TM-score >= 0.5 (or the user's preferred
  threshold).
- USalign route: `usalign_results.jsonl` and a USalign-based candidate list
  with TM-score >= 0.5.
- DEER route: `similarity.csv` from `do_retrieval.py`, sorted by Euclidean
  embedding distance. Lower distance indicates higher structural/functional
  similarity.

## Steps

Choose one route. They are not sequential.

### Step 1: Foldseek Large-Scale Screening

Use `foldseek easy-search` to compare candidate structures against a large
Foldseek database. This route is best when the candidate set is very large.

```bash
foldseek easy-search \
  /path/to/candidate_structures/ \
  /path/to/Foldseek_database/ \
  foldseek_hits.m8 \
  tmp_folder \
  --exhaustive-search 1 \
  -s 7.5 \
  --seq-id-mode 2 \
  --format-output query,target,pident,nident,qcov,qlen,tcov,tlen,evalue,qtmscore,ttmscore,alntmscore,lddt,prob
```

Key output fields:

- `qcov`: query coverage, fraction of the query structure aligned.
- `tcov`: target coverage, fraction of the database target structure aligned.
- `qtmscore`: TM-score normalized by query length.
- `ttmscore`: TM-score normalized by target length.
- `alntmscore`: alignment TM-score.

Purpose:

- Supports fast large-scale screening.
- When this route is used, Foldseek TM-score is the final threshold for the
  shortlist unless the user asks for a strict USalign confirmation in addition.

### Step 2: USalign Fine Structural Alignment

Run USalign pairwise directly on the candidate structures against each template
host enzyme structure. This route is preferred when the candidate set is small
or the user explicitly wants USalign evidence.

```bash
USalign query_structure.pdb \
  -dir2 candidate_structures/ \
  -split 0 \
  candidate_list \
  -outfmt 2 > usalign_tmscores.txt
```

Filtering:

- Keep only candidate proteins with USalign `TM-score >= 0.5`.
- Write `usalign_results.jsonl` for the upstream evidence chain.
- The final candidate list contains only USalign-passing candidates.

### Step 3: DEER Comparison (CPU Only)

Use DEER for structure-based dense enzyme retrieval when the user requests it.
The default execution mode for this project is CPU only. If both the server has a
GPU and the user explicitly confirms they want to use it, they may opt into GPU;
otherwise keep CPU mode.

Repository: [WangJiuming/deer](https://github.com/WangJiuming/deer)

Install the CPU-only environment:

```bash
git clone https://github.com/WangJiuming/deer.git
cd deer
conda env create --name deer --file env/env_cpu.yml
conda activate deer
```

The model checkpoint is downloaded from Hugging Face to `./ckpt` automatically.
If auto-download fails:

```bash
huggingface-cli download cuhkaih/deer --local-dir ./ckpt
```

Prepare the two PDB directories:

- `template_pdb_dir`: template host enzyme structures.
- `database_pdb_dir`: candidate bacterial microbe structures.

Run retrieval with the default CPU-only mode:

```bash
CUDA_VISIBLE_DEVICES="" python do_retrieval.py \
  --template_pdb_dir /path/to/template_pdb/ \
  --database_pdb_dir /path/to/candidate_pdb/ \
  --output_dir /path/to/output
```

Do not pass `--use_fa` in CPU mode; Flash Attention is a GPU-only accelerator and
is not used when `CUDA_VISIBLE_DEVICES=""`.

Optional GPU override (only when the user explicitly asks for it):

```bash
conda activate deer
CUDA_VISIBLE_DEVICES="0" python do_retrieval.py \
  --template_pdb_dir /path/to/template_pdb/ \
  --database_pdb_dir /path/to/candidate_pdb/ \
  --output_dir /path/to/output
```

If the user asks for GPU, first confirm the active DEER environment was installed
with a GPU-enabled YAML (`env/env_gpu.yml` or `env/env_gpu_fa.yml`), not
`env/env_cpu.yml`. Passing `--use_fa` is valid only when the user requests it and
the installed environment supports Flash Attention.

Output:

- `similarity.csv` with columns `template_id`, `bacteria_id`, `distance`.
- Lower `distance` means higher similarity.

Recommended shortlisting:

- Keep candidates whose `distance` is below the user-specified cutoff or within
  the top `--top_k` results. If no threshold is given, report the ranked table
  and let the user choose the cutoff instead of inventing one silently.

Notes:

- CPU-only inference is significantly slower than GPU; run it for a bounded
  candidate set or request approval before launching on a very large set. Keep
  CPU as the default unless the user explicitly chooses GPU.
- DEER distance and USalign/Foldseek TM-score are different scoring scales; do
  not convert them directly into each other.

## Choosing Which Route to Use

- Use Step 1 (Foldseek) for large-scale candidate structure sets where exhaustive
  pairwise USalign would be slow.
- Use Step 2 (USalign) for small-to-medium candidate sets or when USalign is
  already available and expected as the final structural evidence.
- Use Step 3 (DEER) when the user explicitly wants dense enzyme retrieval or the
  workflow calls for DEER evidence.
- If the user asks for multiple routes or explicitly wants one method to
  pre-screen before another, run them in that order. Otherwise, use only one
  route and do not claim the other route was run.
