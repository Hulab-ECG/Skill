#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"
cd "$WORKDIR"
mkdir -p "$CLEANDATA_DIR"
activate_env "$CONDA_RNASEQ"

echo "========================================"
echo "[Step 2] Quality Control (fastp)"
echo "Started: $(date)"
echo "========================================"

FASTP=$(which fastp)
for SampleID in $(control_sample_list); do
    R1="${RAWDATA_DIR}/${SampleID}_1.fastq.gz"
    R2="${RAWDATA_DIR}/${SampleID}_2.fastq.gz"
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "  WARNING: raw reads missing for $SampleID (R1=$R1, R2=$R2), skipping QC"
        continue
    fi
    echo ">>> $(date) Processing: ${SampleID}"
    ${FASTP} -i "$R1" \
        -o "${CLEANDATA_DIR}/${SampleID}.clean_R1.fastq.gz" \
        -I "$R2" \
        -O "${CLEANDATA_DIR}/${SampleID}.clean_R2.fastq.gz" \
        -h "${CLEANDATA_DIR}/${SampleID}.filt.html" \
        -j "${CLEANDATA_DIR}/${SampleID}.filt.json" \
        $FASTP_PARAMS -w $THREADS_QC
    echo ">>> $(date) Done: ${SampleID}"
done

echo "[2.1] $(ls "$CLEANDATA_DIR"/*.clean_R1.fastq.gz | wc -l) samples processed"
conda deactivate
echo "========================================"
echo "[Step 2] Complete: $(date)"
echo "========================================"

# Note: this step only QC's the CONTROL_SAMPLES from config.sh.
