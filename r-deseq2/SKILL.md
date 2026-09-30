---
name: r-deseq2
description: "Run the standardized R DESeq2 differential expression workflow. Use when the user wants differential expression analysis, raw count normalization, volcano/MA plots, DEG tables, candidate gene subsets, or UGT/target gene filtering from a count matrix and sample groups."
---

# R-DEseq2

The R-DEseq2 skill wraps a standardized DESeq2 R script that takes a raw count matrix, sample metadata, and one or more contrasts, then produces normalized counts, differential expression tables, and publication-ready PDF plots.

## 1. Run the workflow

The core script lives at `scripts/run_deseq2.R` in this skill directory and expects R packages `DESeq2` and `pheatmap`.

```bash
Rscript <skill_directory>/scripts/run_deseq2.R \
  --count-matrix count_matrix.csv \
  --meta         sample_meta.tsv \
  --formula      ~group \
  --contrasts    blue,control,blue_vs_control \
  --contrasts    red,control,red_vs_control \
  --contrasts    red,blue,red_vs_blue \
  --outdir       deseq2_results
```

Add `--genes candidate_genes.txt` when candidate-gene subsets are needed.

## 2. Required inputs

Use `--count-matrix` for a CSV raw count matrix with gene IDs in the first column and one sample per column. Values should be integer raw counts.

```csv
gene,CTL_1,CTL_2,CTL_3,B70_1,B70_2,B70_3
MSTRG.1,869,1486,990,2913,2666,2763
```

Use `--meta` for a TSV with `sample` and `group` columns. `sample` should match count-matrix column names. Put the control group first when possible.

```text
sample	group
CTL_1	control
CTL_2	control
B70_1	blue
B70_2	blue
```

Use `--formula` for the DESeq2 design. The default is `~group`.

Use one or more `--contrasts` arguments in `numerator,denominator[,label]` format. If no label is given, the script uses `numerator_vs_denominator`.

## 3. Main parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--count-matrix` | required | Raw count CSV path |
| `--meta` | required | Sample metadata TSV path with `sample`, `group` |
| `--formula` | `~group` | DESeq2 design formula |
| `--contrasts` | required | Compare as `numerator,denominator[,label]`; repeat for multiple comparisons |
| `--outdir` | `deseq2_results` | Output directory |
| `--fdr` | `0.05` | Adjusted p-value threshold |
| `--log2fc` | `1.0` | Absolute log2 fold-change threshold |
| `--min-count` | `10` | Low-count filter: minimum value in a sample |
| `--min-samples` | `2` | Min samples above `--min-count` to keep a gene |
| `--genes` | none | Candidate gene list, one ID per line |
| `--alpha` | same as `--fdr` | Alpha passed to DESeq2 `results()` |
| `--shrink` | `none` | LFC shrinkage: `none`, `normal`, `apeglm`, or `ashr` |

When using `--shrink apeglm`, the R package `apeglm` must be installed.

## 4. Output files

All files are written to `--outdir`.

| File | Content |
| --- | --- |
| `normalized_counts.csv` | DESeq2 normalized count matrix |
| `deseq2_<label>_all.csv` | Full result table for each contrast |
| `deseq2_<label>_cand.csv` | Candidate-gene subset for each contrast |
| `significant_cand_union.csv` | Candidate genes significant in any contrast |
| `de_summary.txt` | Text summary of significant genes per contrast |
| `volcano_<label>.pdf` | Volcano plot per contrast |
| `ma_<label>.pdf` | MA plot per contrast |
| `deseq2_run_meta.csv` | Record of parameters used for the run |

Each result table includes `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, and `significant`.

## 5. Practical notes

- The script aligns metadata samples to the count matrix automatically.
- Repeat `--contrasts` for each comparison; this avoids shell tokenization ambiguity.
- Candidate-gene subset outputs are only written when `--genes` is provided.
- `significant` means the gene passes both the FDR and absolute log2FC thresholds.

Detailed reference notes are in [references/deseq2_standard.md](references/deseq2_standard.md).
