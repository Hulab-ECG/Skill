# R-WGCNA

Standardized R WGCNA co-expression workflow for module detection and candidate
gene ranking against a bait gene.

## Install

You can provide your own R environment. The script needs these R packages:

- `WGCNA`
- `pheatmap`
- `fastcluster`

Create an environment with conda if you do not have your own:

```bash
conda env create -f environment.yml
conda activate r-wgcna
```

Or install with your existing R:

```r
install.packages(c("WGCNA", "pheatmap", "fastcluster"))
```

## Quick run

```bash
Rscript /Users/lily/.codex/skills/r-wgcna/scripts/run_wgcna.R   --expr gene_TPM_matrix.tsv   --meta sample_meta.tsv   --bait MSTRG.7391   --cands candidate_genes.txt   --outdir wgcna_results
```

## Main outputs

- `gene_module_assignments.csv`
- `coexpression_partners.txt`
- `candidate_*_pearson.csv`
- `module_dendrogram.pdf`, `module_trait_heatmap.pdf`, and candidate heatmaps
- `bait_module_hub_genes.csv`

See `SKILL.md` for the full parameter table and output list.
