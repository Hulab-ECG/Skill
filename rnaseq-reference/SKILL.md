---
name: rnaseq-reference
description: Run or deploy the RNAseq_reference pipeline for reference-based RNA-seq quantification with HISAT2 and StringTie (count and TPM matrices).
---

# RNAseq Reference

This skill documents and bundles the `RNAseq_reference` pipeline. The pipeline is a reference-based RNA-seq quantification workflow intended for projects where a reference genome and annotations may exist.

Use this skill when the user asks to:

- run, inspect, deploy, or package the `RNAseq_reference` pipeline;
- run reference-based RNA-seq quantification with HISAT2 + StringTie `-eB`;
- build gene count and gene TPM matrices from paired-end FASTQ data;
- decide whether a StringTie assembly / TransDecoder branch is needed for a reference genome with missing annotations.

Do not use this skill for de novo transcriptome assembly unless requested as mode B of this pipeline.

## Pipeline Overview

The pipeline runs five steps in order:

1. `01_download_rawdata.sh`
   Check `00.rawdata/<SAMPLE>_1.fastq.gz` / `_2.fastq.gz`. If paired FASTQ for a sample is missing, download it from NCBI SRA with `prefetch` + `fasterq-dump`.

2. `02_fastp_qc.sh`
   Run `fastp` for adapter trimming, low-quality filtering, and minimum read length 50. Output goes to `01.cleandata/`.

3. `03_prepare_reference.sh`
   Require a local reference FASTA in `02.ref_genome/` (or `.fna.gz`). If a `genomic.gtf` or `genomic.gff` exists, make sure a GTF is available for StringTie, optionally converting GFF with `gffread -T`, extract reference transcripts/proteins, and build the HISAT2 index.

4. `04_hisat2_align.sh`
   Align clean paired reads with `hisat2 --dta`, sort with `samtools sort`, and index each sample BAM. This step does not create a merged `all_samples.sorted.bam`; StringTie quantification only needs each sample BAM.

5. `05_stringtie_quant.sh`
   Auto-detect annotation mode:
   - Mode A: reference GTF/GFF exists. Run `stringtie <sample>.sorted.bam -p 40 -G ref.gtf -eB` per sample, then `prepDE.py` and `build_tpm_matrix.py`.
   - Mode B: no reference GTF/GFF exists. Run per-sample StringTie assembly, then `stringtie --merge`, produce transcripts with `gffread`, run TransDecoder for peptide/CDS, then quantify all samples against `merged.gtf` with `-eB`.

Both modes finish by creating:
- `04.stringtie/gene_count_matrix.csv`
- `04.stringtie/transcript_count_matrix.csv`
- `04.stringtie/gene_TPM_matrix.tsv`

## Key Design Rules

- Default input is local paired-end FASTQ in `00.rawdata/`. SRA download is a fallback only when files are missing.
- Reference files are user-supplied: edit `config.sh` to point `REF_ASSEMBLY`, `REF_GENOME_NAME`, `REF_GTF_GZ`, and `REF_GFF_GZ` at your own assembly files. Example values such as `GCA_014636745.1_...` are not defaults.
- The full pipeline is bundled under `pipeline/RNAseq_reference/` inside this skill.
- Only one conda environment is used: `rnaseq`. It must include `fastp`, `stringtie`, `hisat2`, `samtools`, `gffread`, `seqkit`, `sra-tools`, `seqtk`, `transdecoder`, and `gffcompare`. See `pipeline/RNAseq_reference/envs/env_rnaseq.yml`.
- Deploy by copying `pipeline/RNAseq_reference/`, editing `config.sh` for `CONDA_BASE` / `WORKDIR` / sample names, then running `bash install_envs.sh` and `bash run_pipeline.sh`.
- Existing outputs are reused: each step skips samples or files that are already present.

## Running the Bundled Pipeline

```bash
cd <skill_directory>/pipeline/RNAseq_reference

bash install_envs.sh
# edit config.sh if needed
bash run_pipeline.sh
```

## Deploying for Another User

If you want to distribute only the pipeline, tar the bundled directory:

```bash
tar --exclude='.DS_Store' -czf RNAseq_reference.tar.gz RNAseq_reference
```

Recipients extract it, edit `config.sh`, then run `install_envs.sh` and `run_pipeline.sh`.
