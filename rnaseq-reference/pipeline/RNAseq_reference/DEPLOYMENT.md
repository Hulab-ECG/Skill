# RNAseq_reference 快速部署指南

## 一句话概览

装好 conda → 创建环境 → 改 `config.sh` 两个变量 → 把参考基因组（和可选注释）放进 `02.ref_genome/`、把下机 FASTQ 放进 `00.rawdata/`（可省略，缺了会自动下载）→ 跑 `run_pipeline.sh`。

本流程只做有参 RNA-seq 定量，不做 UGT 注释、不做 WGCNA。

## 流程模式

- 模式 A：`02.ref_genome/` 里有 `genomic.gtf` 或 `genomic.gff`
  直接用参考注释做 `StringTie -eB -G ref.gtf` 定量。

- 模式 B：`02.ref_genome/` 里既没有 `gtf` 也没有 `gff`
  Step 5 自动切换为 StringTie assembly + merge -> TransDecoder (PEP/CDS) -> StringTie `-eB` 定量。

## 1. 前置要求

- Linux 服务器，已安装 Miniconda（能访问 `conda-forge` / `bioconda`）
- 参考基因组 FASTA 文件，放在 `02.ref_genome/` 下
  - 文件名由你在 `config.sh` 中指定，例如 `GCA_xxxx.1_genomic.fna`
  - 也可以是 `.fna.gz`，Step 3 会自动解压
- 可选的注释文件（有任一即可进入模式 A）：
  - 你自己的参考 GTF 文件，例如 `GCA_xxxx.1_genomic.gtf.gz` 或解压后的 `.gtf`
  - 你自己的参考 GFF 文件，例如 `GCA_xxxx.1_genomic.gff.gz` 或解压后的 `.gff`
- 如果走模式 B，还需要 `transdecoder`，`env_rnaseq.yml` 已包含该软件
- Step 1 默认只检查本地 FASTQ；只有某个样本在 `00.rawdata/` 下缺少 `_1.fastq.gz` / `_2.fastq.gz` 时，才用 `prefetch + fasterq-dump` 自动下载 NCBI SRA

## 2. 创建环境

```bash
bash install_envs.sh
```

会自动创建/更新：

- `rnaseq`：fastp, stringtie, hisat2, samtools, gffread, seqkit, transdecoder 等

如果环境名和你服务器上已有环境冲突，先手动改 `config.sh` 里的 `CONDA_RNASEQ`。

## 3. 唯一要改的配置

编辑 `config.sh`，一般只改：

```bash
CONDA_BASE=""   # 例如你的 conda 安装目录；留空则自动探测
WORKDIR=""      # 例如你的工作目录；留空则用本 pipeline 目录
```

其余可选：线程数、样本清单、参考基因组文件名。

## 4. 放置参考文件

```bash
mkdir -p 02.ref_genome
cp /你的路径/genomic.fna 02.ref_genome/
cp /你的路径/genomic.gtf.gz 02.ref_genome/  # 可选
cp /你的路径/genomic.gff.gz 02.ref_genome/  # 可选
```

## 5. 放置/自动下载原始数据

如果你已经有下机双端 FASTQ：

```bash
mkdir -p 00.rawdata
cp /你的路径/SRRxxx_1.fastq.gz 00.rawdata/
cp /你的路径/SRRxxx_2.fastq.gz 00.rawdata/
```

如果某个样本没给 FASTQ，Step 1 会根据 `CONTROL_SAMPLES` 里的 SRR 编号自动下载。

## 6. 跑流程

```bash
bash run_pipeline.sh
```

或分步跑（每步可断点续跑）：

```bash
bash run_pipeline.sh 2   # fastp QC
bash run_pipeline.sh 3   # 参考基因组 + HISAT2 index（可先准备好 00.rawdata / 01.cleandata）
bash run_pipeline.sh 4   # HISAT2 比对 + 每样本 sorted BAM
bash run_pipeline.sh 5   # StringTie 定量（自动选择 -eB 或 assembly+merge+TransDecoder）
```

如果 rawdata 没放：

```bash
bash run_pipeline.sh 1   # 只检查/下载原始 FASTQ
```

长任务后台跑：

```bash
nohup bash run_pipeline.sh > run.log 2>&1 &
```

## 7. 主要输出

```text
00.rawdata/<SampleID>_1/2.fastq.gz          原始 FASTQ（可自动下载；已存在则跳过）
01.cleandata/<sample>.clean_R1/R2.fastq.gz  QC 后 reads
02.ref_genome/hisat2_index.*               HISAT2 index
02.ref_genome/*.genomic.fna                参考基因组
02.ref_genome/*.genomic.gtf/.gff           参考注释（可选）
03.hisat2/<sample>.sorted.bam              每样本 sorted BAM
04.stringtie/gene_count_matrix.csv         基因 x 样本 counts
04.stringtie/gene_TPM_matrix.tsv           基因 x 样本 TPM
05.annotation/transcripts.fasta            模式 B 转录本
05.annotation/*.transdecoder.pep/cds       模式 B 预测蛋白/CDS
```

## 8. 常见问题

- `conda base not found`：在 `config.sh` 手动写死 `CONDA_BASE`。
- `TransDecoder wrapper not found`：确认 `transdecoder` 已装进 `rnaseq` 环境，或手动设置 `TRANS_DECODER_WRAP`。
- 下载/中断后重跑对应 step 即可，已有输出会自动跳过。
