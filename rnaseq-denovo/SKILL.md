---
name: rnaseq-denovo
description: Run or deploy the RNAseq_denovo pipeline for de novo RNA-seq quantification with Trinity assembly, TransDecoder ORF prediction, and Salmon count/TPM matrices.
---

# RNAseq Denovo

This skill documents and bundles the `RNAseq_denovo` pipeline. It is a de novo RNA-seq quantification workflow for projects without a reference genome or annotation, using paired-end FASTQ.

Use this skill when the user asks to:

- run, inspect, deploy, or package the `RNAseq_denovo` pipeline;
- run de novo RNA-seq quantification with Trinity + TransDecoder + Salmon;
- build gene count and gene TPM matrices from paired-end FASTQ without a reference genome;
- extend or adapt the shared raw-data and fastp steps from the reference pipeline into a de novo workflow.

Do not use this skill for reference-based workflows such as HISAT2 + StringTie; use `rnaseq-reference` instead.

## Pipeline Overview

The pipeline runs five steps in order:

1. `01_download_rawdata.sh`
   Check `00.rawdata/<SAMPLE>_1.fastq.gz` / `_2.fastq.gz`. If paired FASTQ for a sample is missing, download it from NCBI SRA with `prefetch` + `fasterq-dump`.

2. `02_fastp_qc.sh`
   Run `fastp` for adapter trimming, low-quality filtering, and minimum read length 50. Output goes to `01.cleandata/`.

3. `03_trinity_assembly.sh`
   Sort clean reads by name with `seqkit sort -n`, then run Trinity with read normalization. Output is `02.trinity/trinity_out_dir/Trinity.fasta`.

4. `04_transdecoder.sh`
   Run `TransDecoder.LongOrfs` and `TransDecoder.Predict` on Trinity transcripts. Output is `03.transdecoder/Trinity.fasta.{pep,cds,gff3}`.

5. `05_salmon_quant.sh`
   Build a Salmon index, quantify each sample, then run `abundance_estimates_to_matrix.pl` to merge gene-level matrices:
   - `04.salmon/matrix/Trinity.gene.counts.matrix`
   - `04.salmon/matrix/Trinity.gene.TMM.TPM.matrix`

## Key Design Rules

- Default input is local paired-end FASTQ in `00.rawdata/`. SRA download is a fallback only when files are missing.
- The full pipeline is bundled under `pipeline/RNAseq_denovo/` inside this skill.
- Only one conda environment is used: `rnaseq`. It must include `fastp`, `trinity`, `transdecoder`, `salmon`, `hmmer`, `seqkit`, `samtools`, `sra-tools`, and `seqtk`. See `pipeline/RNAseq_denovo/envs/env_rnaseq.yml`.
- `config.sh` controls `CONDA_BASE`, `WORKDIR`, `CONTROL_SAMPLES`, threads, memory, Trinity flags, and Salmon parameters. Common edits are `CONDA_BASE`, `WORKDIR`, and `CONTROL_SAMPLES`.
- Existing outputs are reused: each step skips samples or files that are already present.
- The pipeline is portable: copy the `pipeline/RNAseq_denovo/` directory to the target server, edit `pipeline/RNAseq_denovo/config.sh`, then run `bash pipeline/RNAseq_denovo/install_envs.sh` and `bash pipeline/RNAseq_denovo/run_pipeline.sh`.

## Running the Bundled Pipeline

```bash
cd <skill_directory>/pipeline/RNAseq_denovo

bash install_envs.sh
# edit config.sh if needed
bash run_pipeline.sh
```

## Deploying for Another User

If you want to distribute only the pipeline, tar the bundled directory:

```bash
tar --exclude='.DS_Store' -czf RNAseq_denovo.tar.gz RNAseq_denovo
```

Recipients extract it, edit `config.sh`, then run `install_envs.sh` and `run_pipeline.sh`.
