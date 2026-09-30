#!/bin/bash
#=========================================#
# Step 4: HISAT2 alignment to reference   #
#=========================================#
# Aligns clean paired-end reads to reference genome
# Produces sorted BAM files for each sample.
# StringTie quantification only needs each sample's own sorted BAM.

set -euo pipefail
export ZSH_VERSION=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"
cd "$WORKDIR"

activate_env "$CONDA_RNASEQ"

echo "========================================"
echo "[Step 4] HISAT2 Alignment"
echo "Started: $(date)"
echo "========================================"

HISAT2=$(which hisat2)
SAMTOOLS=$(which samtools)
mkdir -p "$ALIGN_DIR"

# ---- 4.1: Align each sample ----
echo "[4.1] Aligning samples to reference genome..."
for SampleID in $(control_sample_list); do
    R1="${CLEANDATA_DIR}/${SampleID}.clean_R1.fastq.gz"
    R2="${CLEANDATA_DIR}/${SampleID}.clean_R2.fastq.gz"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "  WARNING: clean reads missing for $SampleID (R1=$R1, R2=$R2), skipping"
        continue
    fi

    if [ -f "${ALIGN_DIR}/${SampleID}.sorted.bam" ]; then
        echo "  Skipping $SampleID (already aligned)"
        continue
    fi

    echo "  >>> $SampleID ($(date))"

    $HISAT2 -x "$HISAT2_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS_HISAT2" \
        --dta \
        --summary-file "${ALIGN_DIR}/${SampleID}.hisat2.log" \
        2> "${ALIGN_DIR}/${SampleID}.hisat2.err" \
        | $SAMTOOLS sort -@ "$THREADS_HISAT2" -o "${ALIGN_DIR}/${SampleID}.sorted.bam" -

    $SAMTOOLS index "${ALIGN_DIR}/${SampleID}.sorted.bam"
    echo "  <<< $SampleID done ($(date))"
done

# ---- 4.2: Alignment summary ----
echo "[4.2] Alignment Summary:"
for bam in "${ALIGN_DIR}"/*.sorted.bam; do
    s=$(basename "$bam" .sorted.bam)
    mapped=$($SAMTOOLS flagstat "$bam" | grep "properly paired" | awk '{print $1}' || true)
    total=$($SAMTOOLS flagstat "$bam" | head -1 | awk '{print $1}' || true)
    if [ -n "$total" ] && [ "$total" -gt 0 ] 2>/dev/null; then
        pct=$(python3 -c "print(f'{$mapped/$total*100:.1f}%')" 2>/dev/null || echo "N/A")
    else
        pct="N/A"
    fi
    echo "  ${s}: ${mapped}/${total} properly paired (${pct})"
done

echo "========================================"
echo "[Step 4] Complete: $(date)"
echo "========================================"
conda deactivate
