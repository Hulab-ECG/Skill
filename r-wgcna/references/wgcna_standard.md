# R-WGCNA Standard Workflow Reference

Source workflow: an internal project snapshot. The maintained script is copied into this skill at `scripts/run_wgcna.R`.

## Directory layout of the original standard workflow

```text
wgcna_standard/
  README.md
  run_wgcna.R
  example_input/
    example_TPM_matrix.tsv
    example_sample_meta.tsv
    example_candidate_genes.txt
  example_output/
```

## Input file format details

### Expression matrix

TSV with first column `gene_id` and one sample column per sample. All sample columns must be numeric.

```text
gene_id	SRR12640679	SRR12640678	...
MSTRG.1	12.5	13.1	...
MSTRG.2	510.2	520.4	...
```

### Sample metadata

TSV with at least two columns: `sample` and `trait`.

```text
sample	trait
SRR12640679	leaf
SRR12639761	root
SRR12649935	stem
```

The script keeps only metadata samples present in the expression matrix and binaries traits for module-trait correlation.

### Candidate gene list

Optional plain text file with one gene ID per line, often an HMMER PF00201-predicted UGT list.

```text
MSTRG.1
MSTRG.7391
```

## Example run

```bash
Rscript <skill_directory>/scripts/run_wgcna.R \
  --expr      example_input/example_TPM_matrix.tsv \
  --meta      example_input/example_sample_meta.tsv \
  --bait      MSTRG.7391 \
  --cands     example_input/example_candidate_genes.txt \
  --outdir    example_output
```

Full parameter example:

```bash
Rscript <skill_directory>/scripts/run_wgcna.R \
  --expr        example_input/example_TPM_matrix.tsv \
  --meta        example_input/example_sample_meta.tsv \
  --bait        MSTRG.7391 \
  --cands       example_input/example_candidate_genes.txt \
  --outdir      example_output \
  --network     unsigned \
  --power       6 \
  --min-size    30 \
  --merge       0.25 \
  --max-block   30000 \
  --threshold   1 \
  --min-samples 5 \
  --nthreads    4
```

## Parameter table

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--expr` | required | Gene expression matrix path |
| `--meta` | required | Sample metadata path |
| `--bait` | required | Bait gene ID |
| `--cands` | optional | Candidate gene list |
| `--outdir` | `wgcna_results` | Output directory |
| `--network` | `unsigned` | `unsigned`, `signed`, or `signed hybrid` |
| `--power` | auto | Fixed soft threshold; omit for auto selection |
| `--min-size` | `30` | Minimum module size |
| `--merge` | `0.25` | Module merge cut height |
| `--max-block` | `30000` | Max block size for `blockwiseModules` |
| `--threshold` | `1` | Gene low-expression filter threshold (TPM) |
| `--min-samples` | auto | Minimum number of samples above threshold |
| `--nthreads` | `10` | WGCNA threads |
| `--save-tom` | `0` | Save TOM, `1` yes / `0` no |

## Output files

```text
wgcna_summary.csv                  运行摘要
sample_clustering.pdf              样本聚类树
soft_threshold.pdf                 软阈值选择图
module_dendrogram.pdf              模块聚类图
module_trait_heatmap.pdf           模块-性状热图
module_trait_correlation.csv       模块-性状相关和 p 值
gene_module_assignments.csv        所有基因的模块归属
coexpression_partners.txt          bait 模块全部共表达基因
ugt_in_bait_module.txt             bait 模块中的 UGT 基因
candidate_ugt_pearson.csv          候选 UGT 与 bait 的 Pearson 排序
ugt_to_module_map.csv              UGT 到模块的对应关系
ugt_tpm_matrix.csv                 UGT 基因的表达子矩阵
ugt_bait_module_heatmap.pdf        bait 模块候选 UGT 表达热图
bait_module_hub_genes.csv          bait 模块按 kME 排序的 hub 基因
```

## R dependencies

```r
WGCNA
pheatmap
fastcluster
```

## Known caveats

- The example matrix is ~2000 simulated genes with leaf/root/stem patterns and usually produces multiple real modules.
- Small matrices may produce only the `grey` module.
- `mergeCloseModules` may print `Color levels are empty` on tiny or low-signal datasets; this is expected for datasets too small to form co-expression modules.
