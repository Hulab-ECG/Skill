---
name: r-wgcna
description: "Run the standardized R WGCNA co-expression workflow for candidate gene mining. Use when the user wants module detection, bait-module coexpression partners, hub genes, module-trait correlations, or UGT/other candidate ranking from a gene expression matrix."
---

# R-WGCNA

The R-WGCNA skill wraps a standardized WGCNA R script that takes a gene expression matrix, sample metadata, and a bait gene, then produces fixed module/co-expression outputs. This is designed for workflows like UGT candidate mining, but it works for any gene set where WGCNA co-expression analysis is the goal.

## 1. Run the workflow

The core script lives at `scripts/run_wgcna.R` in this skill directory and expects R packages `WGCNA`, `pheatmap`, and `fastcluster`.

```bash
Rscript <skill_directory>/scripts/run_wgcna.R \
  --expr  gene_TPM_matrix.tsv \
  --meta  sample_meta.tsv \
  --bait  MSTRG.7391 \
  --cands ugt_candidates.txt \
  --outdir wgcna_results
```

If only a network is needed and no candidate list exists, omit `--cands`.

## 2. Required inputs

Use `--expr` for a TSV expression matrix with the first column named `gene_id` and one column per sample, all sample columns numeric.

```text
gene_id	SRR12640679	SRR12640678
MSTRG.1	12.5	13.1
MSTRG.2	510.2	520.4
```

Use `--meta` for a TSV with at least `sample` and `trait` columns. The script intersects metadata samples with expression-matrix columns and binaries the trait for module-trait correlation.

Use `--bait` for the bait gene ID, which must be present in the expression matrix.

Use `--cands` for an optional plain-text list of candidate gene IDs, one per line. They are ranked against the bait by Pearson correlation and mapped to WGCNA modules.

## 3. Main parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `wgcna_results` | Output directory |
| `--network` | `unsigned` | Network type: `unsigned`, `signed`, or `signed hybrid` |
| `--power` | auto | Fixed soft threshold; omit for automatic selection |
| `--min-size` | `30` | Minimum module size |
| `--merge` | `0.25` | Module merge cut height |
| `--max-block` | `30000` | Max genes per block in `blockwiseModules` |
| `--threshold` | `1` | Expression filter threshold, typically TPM |
| `--min-samples` | auto | Min samples above threshold required to keep a gene |
| `--nthreads` | `10` | WGCNA thread count |
| `--save-tom` | `0` | Set `1` to save TOM output |

## 4. Output files

All results are written to `--outdir`.

| File | Content |
| --- | --- |
| `wgcna_summary.csv` | Run settings and key summary metrics |
| `sample_clustering.pdf` | Sample hierarchical clustering and trait heatmap |
| `soft_threshold.pdf` | Soft-threshold selection plots |
| `module_dendrogram.pdf` | Gene dendrogram and module colors |
| `module_trait_heatmap.pdf` | Module-trait heatmap |
| `module_trait_correlation.csv` | Module-trait correlation and p-value table |
| `gene_module_assignments.csv` | Gene-to-module assignments for all genes |
| `coexpression_partners.txt` | All genes in the bait gene's module |
| `ugt_in_bait_module.txt` | Candidate genes found in the bait module |
| `candidate_ugt_pearson.csv` | Candidate genes ranked by Pearson correlation to the bait |
| `ugt_to_module_map.csv` | Candidate gene-to-module mapping |
| `ugt_tpm_matrix.csv` | Candidate gene expression submatrix |
| `ugt_bait_module_heatmap.pdf` | Candidate gene expression heatmap in the bait module |
| `bait_module_hub_genes.csv` | Bait-module hub genes ranked by kME |

## 5. Practical notes

- The script filters low-expression genes using `--threshold` and `--min-samples` before WGCNA.
- For small test matrices, WGCNA often detects only the `grey` module. That is usually a data-size limitation, not a script error.
- `mergeCloseModules` can print messages such as `Color levels are empty` when the network has no real modules.
- Candidate names on the PDF and CSV outputs keep the generic "UGT" wording from the source script, but the logic applies to any candidate gene list.

Detailed reference notes are in [references/wgcna_standard.md](references/wgcna_standard.md).
