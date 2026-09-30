# R-DEseq2 Standard Workflow Reference

Source workflow: an internal project snapshot. The maintained script is copied into this skill at `scripts/run_deseq2.R`.

## Directory layout of the original standard workflow

```text
deseq2_standard/
  README.md
  run_deseq2.R
  example_input/
    example_count_matrix.csv
    example_sample_meta.tsv
    example_candidate_genes.txt
  example_output/
    de_summary.txt
    volcano_blue_vs_control.pdf
    volcano_red_vs_control.pdf
    volcano_red_vs_blue.pdf
```

## Input file format details

### Raw count matrix

CSV with first column as gene ID and one sample column per sample. Values should be integer raw counts.

```csv
gene,CTL_1,CTL_2,CTL_3,B70_1,B70_2,B70_3,R70_1,R70_2,R70_3
MSTRG.1,869,1486,990,2913,2666,2763,942,868,885
```

The first column name is flexible; the script treats the first column as the gene ID column.

### Sample metadata

TSV with exactly the needed columns `sample` and `group`.

```text
sample	group
CTL_1	control
CTL_2	control
B70_1	blue
B70_2	blue
R70_1	red
R70_2	red
```

`sample` must match count-matrix column names. Extra metadata samples are aligned away automatically.

### Candidate gene list

Optional plain text file with one gene ID per line, often UGT or other target genes.

```text
MSTRG.126
MSTRG.1024
MSTRG.2901
```

When this input is not provided, candidate-specific result files are not written.

## Example run

```bash
Rscript <skill_directory>/scripts/run_deseq2.R \
  --count-matrix example_input/example_count_matrix.csv \
  --meta         example_input/example_sample_meta.tsv \
  --formula      ~group \
  --contrasts    blue,control,blue_vs_control \
  --contrasts    red,control,red_vs_control \
  --contrasts    red,blue,red_vs_blue \
  --genes        example_input/example_candidate_genes.txt \
  --outdir       example_output
```

## Parameter table

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--count-matrix` | required | Raw count CSV path |
| `--meta` | required | Sample metadata TSV path with `sample`, `group` |
| `--formula` | `~group` | DESeq2 design formula; multi-factor can be `~batch+group` |
| `--contrasts` | required | One or more comparisons in `numerator,denominator[,label]` format |
| `--outdir` | `deseq2_results` | Output directory |
| `--fdr` | `0.05` | Adjusted p-value threshold |
| `--log2fc` | `1.0` | Absolute log2 fold-change threshold |
| `--min-count` | `10` | Low-count filter threshold |
| `--min-samples` | `2` | Min samples above threshold to keep a gene |
| `--genes` | none | Candidate gene list path |
| `--alpha` | same as `--fdr` | Alpha passed to DESeq2 `results()` |
| `--shrink` | `none` | LFC shrinkage: `none`, `normal`, `apeglm`, or `ashr` |

## Contrast syntax

Each comparison is written as one `--contrasts` argument:

```bash
--contrasts blue,control,blue_vs_control \
--contrasts red,control,red_vs_control
```

Without an explicit label, the default label is `numerator_vs_denominator`.

The last variable in `--formula` is treated as the comparison factor. With `--formula ~group`, the numerator and denominator must both be levels in `group`.

Using multiple `--contrasts` arguments avoids shell splitting problems.

## Output files

```text
normalized_counts.csv           DESeq2 归一化 count
deseq2_<label>_all.csv          每个比较的全部基因结果
deseq2_<label>_cand.csv          每个比较的候选基因子集
significant_cand_union.csv       至少一个比较显著的候选基因并集
de_summary.txt                   文本结果摘要
volcano_<label>.pdf              每个比较的火山图（矢量图）
ma_<label>.pdf                   每个比较的 MA 图（矢量图）
deseq2_run_meta.csv              本次运行的参数记录
```

Result tables contain `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, and `significant`.

## R dependencies

```r
DESeq2
pheatmap
```

Optional dependency:

```r
apeglm
```

`apeglm` is used only when `--shrink apeglm` is specified.

## Relationship to the original pipeline script

The original hardcoded script was `09_deseq2_light.R` with fixed control/blue/red groups, leaf naming, and UGT naming. This standardized wrapper keeps the core DESeq2 logic but makes grouping, comparison design, and output names generic.
