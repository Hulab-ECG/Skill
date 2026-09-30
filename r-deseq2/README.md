# R-DEseq2

Standardized R DESeq2 differential expression workflow. Use it when you have a
matrix of raw read counts plus sample groups and need normalized counts, DEG
tables, volcano/MA plots, candidate-gene subsets, or target filtering.

## Install

You can provide your own R environment. The script needs these R packages:

- `DESeq2`
- `pheatmap`
- optionally `apeglm` when `--shrink apeglm` is requested

Create an environment with conda if you do not have your own:

```bash
conda env create -f environment.yml
conda activate r-deseq2
```

Or install with your existing R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
install.packages("pheatmap")
```

## Quick run

```bash
Rscript /Users/lily/.codex/skills/r-deseq2/scripts/run_deseq2.R   --count-matrix count_matrix.csv   --meta sample_meta.tsv   --formula ~group   --contrasts blue,control,blue_vs_control   --contrasts red,control,red_vs_control   --outdir deseq2_results
```

## Main outputs

- `normalized_counts.csv`
- `deseq2_<label>_all.csv` and `deseq2_<label>_cand.csv`
- `volcano_<label>.pdf` and `ma_<label>.pdf`
- `de_summary.txt` and `deseq2_run_meta.csv`

See `SKILL.md` for the full parameter table and input conventions.
