#!/bin/bash
#=========================================#
# Step 3: Trinity de novo assembly        #
#=========================================#
# Sorts reads by name to fix R1/R2 pairing, then runs Trinity
# with read normalization and THREADS_TRINITY CPU cores.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"

cd "$WORKDIR"
export LC_ALL=C

activate_env "$CONDA_RNASEQ" || exit 1

echo "========================================"
echo "[Step 3] Trinity De Novo Assembly"
echo "Started: $(date)"
echo "========================================"

# ---- 3.1: Sort reads by name ----
mkdir -p "$SORTED_DIR"
SEQKIT=$(which seqkit)

echo "[3.1] Sorting reads by name (seqkit sort -n)..."
for f in "$CLEANDATA_DIR"/*.clean_R1.fastq.gz; do
    base=$(basename "$f" .clean_R1.fastq.gz)
    echo "  Sorting $base..."

    if [ ! -f "$SORTED_DIR/${base}.clean_R1.fastq.gz" ]; then
        $SEQKIT sort -n -j "$THREADS_SORT" \
            -o "$SORTED_DIR/${base}.clean_R1.fastq.gz" \
            "$CLEANDATA_DIR/${base}.clean_R1.fastq.gz" &
        $SEQKIT sort -n -j "$THREADS_SORT" \
            -o "$SORTED_DIR/${base}.clean_R2.fastq.gz" \
            "$CLEANDATA_DIR/${base}.clean_R2.fastq.gz" &
        wait
    else
        echo "  Skipping $base (already sorted)"
    fi
done
echo "[3.1] Sorting done: $(date)"

# ---- 3.2: Build comma-separated file lists ----
echo "[3.2] Building file lists..."
ls "$SORTED_DIR"/*.clean_R1.fastq.gz 2>/dev/null | sort > "$SORTED_DIR/left_list_sorted.txt" || true
ls "$SORTED_DIR"/*.clean_R2.fastq.gz 2>/dev/null | sort > "$SORTED_DIR/right_list_sorted.txt" || true
LEFT_FILES=$(paste -sd, "$SORTED_DIR/left_list_sorted.txt")
RIGHT_FILES=$(paste -sd, "$SORTED_DIR/right_list_sorted.txt")

if [ -z "$LEFT_FILES" ] || [ -z "$RIGHT_FILES" ]; then
    echo "  ERROR: no sorted paired FASTQ files found in $SORTED_DIR" >&2
    exit 1
fi

# ---- 3.3: Run Trinity ----
echo "[3.3] Running Trinity..."
echo "  CPU: $THREADS_TRINITY"
echo "  Memory: $MAX_MEMORY"
echo "  Normalization: enabled"
echo "  Started: $(date)"

mkdir -p "$(dirname "$TRINITY_DIR")"
Trinity --seqType fq \
        --max_memory "$MAX_MEMORY" \
        --left "$LEFT_FILES" \
        --right "$RIGHT_FILES" \
        --CPU "$THREADS_TRINITY" \
        --output "$TRINITY_DIR" \
        $TRINITY_NORMALIZE \
        $TRINITY_EXTRA

echo "[3.4] Assembly Summary:"
echo "  Transcript count: $(grep -c '>' "$TRINITY_FA" 2>/dev/null || echo 0)"
ls -lh "$TRINITY_FA"

echo "========================================"
echo "[Step 3] Complete: $(date)"
echo "========================================"
conda deactivate || true
