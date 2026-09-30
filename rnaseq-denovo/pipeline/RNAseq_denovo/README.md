# RNAseq_denovo Pipeline

无参 RNA-seq 定量流程：从原始 FASTQ 到 Trinity 转录本、TransDecoder 蛋白/CDS、Salmon 定量及基因表达矩阵。

前两步与 `RNAseq_reference` 共用：
```text
01_download_rawdata.sh
02_fastp_qc.sh
```

后面参考 centella de novo portable 的组装和定量流程：
```text
03_trinity_assembly.sh
04_transdecoder.sh
05_salmon_quant.sh
```

## 流程

```text
输入数据检查/自动下载 (00.rawdata)
  -> fastp QC (01.cleandata)
  -> 按 sample name 排序 reads (01.cleandata.sorted)
  -> Trinity de novo assembly (02.trinity)
  -> TransDecoder ORF prediction (03.transdecoder)
  -> Salmon index + quant (04.salmon)
  -> abundance_estimates_to_matrix.pl 合并 count/TPM 矩阵
```

## 目录

```text
config.sh                                        # 配置：conda、workdir、样本、线程
run_pipeline.sh                                  # 主运行脚本
install_envs.sh                                  # 创建 conda 环境
envs/
  env_rnaseq.yml                                 # fastp, trinity, transdecoder, salmon, seqkit, sra-tools, samtools, hmmer, seqtk
scripts/
  01_download_rawdata.sh                         # 检查 00.rawdata；缺失时自动下载 SRA
  02_fastp_qc.sh                                 # fastp adapter/quality trim
  03_trinity_assembly.sh                         # 排序 reads + Trinity 组装
  04_transdecoder.sh                             # TransDecoder LongOrfs + Predict
  05_salmon_quant.sh                             # salmon index/quant + 表达矩阵
```

## 快速部署

```bash
tar -xzf RNAseq_denovo.tar.gz
cd RNAseq_denovo

# 1. 创建环境（首次需要网络）
bash install_envs.sh

# 2. 按需编辑 config.sh
vi config.sh

# 3. 跑全套
bash run_pipeline.sh

# 或分步跑
bash run_pipeline.sh 1   # 下载/检查原始 FASTQ
bash run_pipeline.sh 2   # fastp QC
bash run_pipeline.sh 3   # Trinity 组装
bash run_pipeline.sh 4   # TransDecoder
bash run_pipeline.sh 5   # Salmon 定量和表达矩阵
```

## 关键输出

```text
00.rawdata/<SampleID>_1/2.fastq.gz          原始双端 FASTQ
01.cleandata/<SampleID>.clean_R1/R2.fastq.gz  QC 后 reads
02.trinity/trinity_out_dir/Trinity.fasta      Trinity 转录本
03.transdecoder/Trinity.fasta.pep              蛋白序列
03.transdecoder/Trinity.fasta.cds              CDS 序列
03.transdecoder/Trinity.fasta.gff3             ORF 注释
04.salmon/matrix/Trinity.gene.counts.matrix     基因 raw counts
04.salmon/matrix/Trinity.gene.TMM.TPM.matrix    基因 TMM-TPM，适合 WGCNA
```

## 说明

- Step 1 默认不下载数据：如果你已经在 `00.rawdata/` 下放好了双端 FASTQ，它会直接跳过下载。
- Step 3 会把 clean reads 按 read name 排序，保证 Trinity 能正确识别 R1/R2 配对。
- Step 5 使用 Salmon 对 trinity 转录本定量，并通过 `abundance_estimates_to_matrix.pl` 生成下游分析常用的 counts/TPM 矩阵。
