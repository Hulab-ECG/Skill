#!/bin/bash
#=========================================#
# Step 1: Check/download raw paired-end   #
#         RNA-seq data                    #
#=========================================#
# Default: user already places paired-end FASTQ files in 00.rawdata/:
#   ${SampleID}_1.fastq.gz
#   ${SampleID}_2.fastq.gz
# If a sample is missing, this script downloads it from NCBI SRA using
# prefetch + fasterq-dump (sra-tools in the rnaseq conda env).

set -euo pipefail
export ZSH_VERSION=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"

cd "$WORKDIR"
mkdir -p "$RAWDATA_DIR"

activate_env "$CONDA_RNASEQ"

echo "========================================"
echo "[Step 1] Raw FASTQ check / SRA download"
echo "Project: $BIOPROJECT"
echo "Started: $(date)"
echo "========================================"

# Build the list of samples that do not yet have both paired FASTQ files.
: > "$RAWDATA_DIR/missing_list.txt"
for SampleID in $(control_sample_list); do
    R1="${RAWDATA_DIR}/${SampleID}_1.fastq.gz"
    R2="${RAWDATA_DIR}/${SampleID}_2.fastq.gz"
    if [ -s "$R1" ] && [ -s "$R2" ]; then
        echo "  OK  $SampleID (local FASTQ exists)"
    else
        echo "$SampleID" >> "$RAWDATA_DIR/missing_list.txt"
        echo "  MISSING  $SampleID"
    fi
done

if [ ! -s "$RAWDATA_DIR/missing_list.txt" ]; then
    echo "[1.1] All control samples already have paired-end FASTQ files."
    echo "      No download is needed."
    echo "========================================"
    echo "[Step 1] Complete: $(date)"
    echo "========================================"
    conda deactivate
    exit 0
fi

echo "[1.1] Samples to fetch from NCBI SRA: $(wc -l < "$RAWDATA_DIR/missing_list.txt")"
echo "      One SRR accession per line in $RAWDATA_DIR/missing_list.txt"

cd "$RAWDATA_DIR"

# ---- 1.2: prefetch missing SRA records ----
echo "[1.2] Downloading SRA files with prefetch..."
for srr in $(cat missing_list.txt); do
    if [ -f "${srr}/${srr}.sra" ] || [ -f "${srr}.sra" ] || [ -f "${srr}_1.fastq.gz" ]; then
        echo "  Skipping download $srr (SRA/FASTQ already exists locally)"
    else
        echo "  Downloading $srr..."
        prefetch "$srr" -o "${srr}/${srr}.sra" --max-size 50G
    fi
done

# ---- 1.3: fasterq-dump + gzip ----
echo "[1.3] Converting SRA to FASTQ with fasterq-dump..."
for srr in $(cat missing_list.txt); do
    if [ -f "${srr}_1.fastq.gz" ]; then
        echo "  Skipping conversion $srr (FASTQ already exists)"
        continue
    fi
    echo "  Converting $srr..."
    fasterq-dump "$srr" --split-3 --threads "$THREADS_DOWNLOAD" --outdir .
    if [ -f "${srr}_1.fastq" ]; then
        gzip -f "${srr}_1.fastq"
    fi
    if [ -f "${srr}_2.fastq" ]; then
        gzip -f "${srr}_2.fastq"
    fi
done

echo "[1.4] Download Summary:"
ls -lh *_1.fastq.gz *_2.fastq.gz 2>/dev/null || true

conda deactivate

echo "========================================"
echo "[Step 1] Complete: $(date)"
echo "========================================"
