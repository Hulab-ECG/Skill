# RNAseq_denovo 快速部署指南

## 一句话概览

装好 conda → 创建 `rnaseq` 环境 → 改 `config.sh` → 把下机双端 FASTQ 放进 `00.rawdata/`（可省略，缺了自动下载）→ 跑 `run_pipeline.sh`。

本流程是无参 RNA-seq：Trinity 组装 + TransDecoder 注释 + Salmon 定量。

## 1. 前置要求

- Linux 服务器，已安装 Miniconda（能访问 `conda-forge` / `bioconda`）
- 双端下机 FASTQ，或提供可以映射到 `CONTROL_SAMPLES` 的 SRR 编号

## 2. 创建环境

```bash
bash install_envs.sh
```

会自动创建/更新一个环境：

- `rnaseq`：fastp, trinity, transdecoder, salmon, hmmer, seqkit, samtools, sra-tools, seqtk

如果环境名和你服务器上已有环境冲突，先手动改 `config.sh` 里的 `CONDA_RNASEQ`。

## 3. 唯一要改的配置

编辑 `config.sh`，一般只改：

```bash
CONDA_BASE=""   # 例如你的 conda 安装目录；留空则自动探测
WORKDIR=""      # 例如你的工作目录；留空则用本 pipeline 目录
```

其余可选：线程数、内存、样本清单、Salmon k-mer。

## 4. 放置原始数据

```bash
mkdir -p 00.rawdata
cp /你的路径/SRRxxx_1.fastq.gz 00.rawdata/
cp /你的路径/SRRxxx_2.fastq.gz 00.rawdata/
```

如果某个样本没给 FASTQ，Step 1 会根据 `CONTROL_SAMPLES` 里的 SRR 编号自动下载。

## 5. 跑流程

```bash
bash run_pipeline.sh
```

或分步跑：

```bash
bash run_pipeline.sh 2   # fastp QC
bash run_pipeline.sh 3   # Trinity 组装（内存和耗时高）
bash run_pipeline.sh 4   # TransDecoder
bash run_pipeline.sh 5   # Salmon 定量 + 表达矩阵
```

## 6. 主要输出

```text
02.trinity/trinity_out_dir/Trinity.fasta         Trinity 转录本
03.transdecoder/Trinity.fasta.pep                蛋白序列
03.transdecoder/Trinity.fasta.cds                CDS 序列
04.salmon/matrix/Trinity.gene.counts.matrix      基因 counts
04.salmon/matrix/Trinity.gene.TMM.TPM.matrix     基因 TMM-TPM 矩阵
```

## 7. 常见问题

- `conda base not found`：在 `config.sh` 手动写死 `CONDA_BASE`。
- `Trinity not found`：确认 `trinity` 已装进 `rnaseq` 环境。
- `abundance_estimates_to_matrix.pl not found`：确认 Trinity 环境的 `PATH` 正常，通常位于 `rnaseq/opt/trinity/.../util/` 或通过 conda 自动暴露。
- 下载/中断后重跑对应 step 即可，已有输出会自动跳过。
