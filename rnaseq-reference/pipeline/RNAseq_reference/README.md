# RNAseq_reference Pipeline

精简版有参 RNA-seq 定量流程：保留原始数据检查/自动下载、QC、参考基因组准备、
HISAT2 比对、StringTie 定量这几个核心环节，不做 UGT 注释、不做 WGCNA、不做候选基因提取。
目标输出是基因水平的 **count 矩阵**和 **TPM 矩阵**，供后续差异表达或共表达分析使用。

默认认为你已经把下机后的双端 FASTQ 文件放在 `00.rawdata/` 下；如果某个样本缺失，Step 1 会根据 `CONTROL_SAMPLES` 里的 SRA 编号自动从 NCBI 下载。

## 流程

```
输入数据检查/自动下载 (00.rawdata)
  -> fastp QC (01.cleandata)
  -> 参考基因组（自动识别 GTF/GFF 是否存在）+ HISAT2 index (02.ref_genome)
  -> HISAT2 --dta 比对 + samtools sort (03.hisat2)
  -> StringTie 定量 (04.stringtie)
       模式 A：有 GTF/GFF -> StringTie -eB -G ref.gtf
       模式 B：无 GTF/GFF -> StringTie assembly + merge + TransDecoder (PEP/CDS)
  -> prepDE.py 合并成 count matrix + build_tpm_matrix.py 合并成 TPM matrix
```

## 目录

```
config.sh                                        # 配置：conda、workdir、样本、参考基因组、线程
run_pipeline.sh                                  # 主运行脚本
install_envs.sh                                  # 创建 conda 环境
envs/
  env_rnaseq.yml                                 # fastp, stringtie, hisat2, gffread, seqkit, samtools, transdecoder
scripts/
  01_download_rawdata.sh                         # 检查 00.rawdata；缺失时用 prefetch+fasterq-dump 下载
  02_fastp_qc.sh                                 # fastp adapter/quality trim
  03_prepare_reference.sh                        # genome + (GTF/GFF, 可选) + HISAT2 index
  04_hisat2_align.sh                             # HISAT2 -> 每样本 sorted BAM
  05_stringtie_quant.sh                          # StringTie -eB 或 assembly+merge+TransDecoder -> count/TPM matrix
  prepDE.py                                      # StringTie -> count matrix helper
  build_tpm_matrix.py                            # gene_abund.tab -> TPM matrix helper
  reference_prot_gene_map.py                     # GFF/GTF ID 映射 helper
  build_prot_gene_map.py                         # merged GTF transcript->gene helper
```

## 快速部署

```bash
tar -xzf RNAseq_reference.tar.gz
cd RNAseq_reference

# 1. 创建环境（首次）
bash install_envs.sh

# 2. 如果 conda 或工作目录特殊，编辑 config.sh
vi config.sh

# 3. 跑全套
bash run_pipeline.sh

# 或分步跑
bash run_pipeline.sh 1   # 只检查/下载原始 FASTQ
bash run_pipeline.sh 3   # 只准备参考基因组 + GTF/GFF + HISAT2 index
bash run_pipeline.sh 5   # 只做 StringTie 定量（自动选择 -eB 或 assembly+merge+TransDecoder）
```

## 样本设置

在 `config.sh` 中把 `CONTROL_SAMPLES` 填成你自己的样本 SRA 编号，每行一个，
例如：

```text
SRR12640679
SRR12640678
```

默认是空字符串；示例中的 `PRJNA642665` / `SRR...` 只是说明格式。

## 关键输出

```text
00.rawdata/<SampleID>_1/2.fastq.gz  原始双端 FASTQ（可自动下载）
04.stringtie/gene_count_matrix.csv    基因 x 样本 raw counts
04.stringtie/gene_TPM_matrix.tsv      基因 x 样本 TPM
03.hisat2/<sample>.sorted.bam         各样本 sorted BAM
02.ref_genome/hisat2_index.*          HISAT2 index
02.ref_genome/*.genomic.fna           reference genome
02.ref_genome/*.genomic.gtf/.gff      reference annotations（可选）
05.annotation/transcripts.fasta      模式 B 生成的转录本序列
05.annotation/*.transdecoder.pep/cds 模式 B 生成的蛋白/CDS
```

## 说明

- Step 1 默认不下载数据：如果你已经在 `00.rawdata/` 下放好了所有 `_1.fastq.gz` / `_2.fastq.gz`，它会直接跳过下载。只有某个样本缺文件时，才用 `prefetch + fasterq-dump` 自动下载对应 SRA。
- Step 3 会先检查 `02.ref_genome/` 下是否有本地参考文件：
  - 你自己的参考基因组文件，例如 `GCA_xxxx.1_genomic.fna`（或 `.fna.gz`）
  - 你自己的参考 GTF，例如 `GCA_xxxx.1_genomic.gtf`（或 `.gtf.gz`）
  - 你自己的参考 GFF，例如 `GCA_xxxx.1_genomic.gff`（或 `.gff.gz`）
  - 如果 GTF/GFF 有任一存在，Step 3 会处理注释并用 `gffread` 提取参考转录本/蛋白；如果只有 GFF，会自动用 `gffread -T` 转成 StringTie 需要的 GTF。
  - 如果 GTF/GFF 都没有，Step 3 只构建普通 HISAT2 index，不报错；Step 5 会自动切换成 StringTie assembly + merge + TransDecoder 模式。
- Step 4 每个样本生成 `*.sorted.bam`，StringTie 定量只用每个样本自己的 sorted BAM。
- Step 5 自动判断模式：
  - 有参考 GTF：对每个样本跑 `StringTie -eB -G ref.gtf`，然后 `prepDE.py` 合并 count matrix、`build_tpm_matrix.py` 合并 TPM matrix。
  - 无参考 GTF：先对每个样本跑 StringTie 组装，再 `stringtie --merge`，用 `gffread` 提转录本，再用 TransDecoder 预测蛋白/CDS，然后基于 `merged.gtf` 跑 `StringTie -eB` 定量。
