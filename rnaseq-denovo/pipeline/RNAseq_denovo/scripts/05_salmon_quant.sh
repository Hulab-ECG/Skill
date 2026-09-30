#!/bin/bash
#=========================================#
# Step 5: Salmon quantification           #
#=========================================#
# Indexes Trinity transcripts and quantifies expression for all samples
# using clean reads. Builds gene-level count and TPM matrices via
# Trinity's abundance_estimates_to_matrix.pl.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"

cd "$WORKDIR"

activate_env "$CONDA_RNASEQ" || exit 1

echo "========================================"
echo "[Step 5] Salmon Quantification"
echo "Started: $(date)"
echo "========================================"

echo "[5.0] Waiting for TransDecoder..."
while [ ! -f "$PEP_FILE" ]; do
    sleep 60
done
sleep 10
echo "[5.0] TransDecoder done: $(date)"

# ---- 5.1: Build Salmon index ----
echo "[5.1] Building Salmon index (k=$SALMON_KMER)..."
mkdir -p "$SALMON_DIR"
salmon index -t "$TRINITY_FA" \
    -p "$THREADS_SALMON" \
    -k "$SALMON_KMER" \
    --type "$SALMON_TYPE" \
    -i "$SALMON_DIR/index"
echo "  Index done: $(date)"

# ---- 5.2: Quantify each sample ----
echo "[5.2] Quantifying samples..."
for filename in "${RAWDATA_DIR}"/*_1.fastq.gz; do
    SampleID=$(basename "$filename" _1.fastq.gz)
    if [ ! -f "${CLEANDATA_DIR}/${SampleID}.clean_R1.fastq.gz" ] || [ ! -f "${CLEANDATA_DIR}/${SampleID}.clean_R2.fastq.gz" ]; then
        echo "  WARNING: clean reads missing for $SampleID, skipping"
        continue
    fi
    if [ -f "${SALMON_DIR}/${SampleID}.quant/quant.sf" ]; then
        echo "  Skipping $SampleID (already quantified)"
        continue
    fi
    echo "  >>> $SampleID ($(date))"

    salmon quant -i "$SALMON_DIR/index" \
        -l A \
        -p "$THREADS_SALMON" \
        --validateMappings \
        -1 "${CLEANDATA_DIR}/${SampleID}.clean_R1.fastq.gz" \
        -2 "${CLEANDATA_DIR}/${SampleID}.clean_R2.fastq.gz" \
        -o "${SALMON_DIR}/${SampleID}.quant"

    echo "  <<< $SampleID done ($(date))"
done

# ---- 5.3: Build gene-level expression matrices ----
echo "[5.3] Building gene-level matrices via abundance_estimates_to_matrix.pl..."
mkdir -p "$MATRIX_DIR"

if [ ! -f "$GENE_TRANS_MAP" ]; then
    echo "ERROR: gene_trans_map not found: $GENE_TRANS_MAP" >&2
    echo "Trinity assembly may be incomplete." >&2
    exit 934
fi

abundance_estimates_to_matrix.pl \
    --est_method salmon \
    --out_prefix "$MATRIX_DIR/Trinity.gene" \
    --name_sample_by_basedir \
    --cross_sample_norm TMM \
    --gene_trans_map "$GENE_TRANS_MAP" \
    "$SALMON_DIR"/*.quant/quant.sf

echo "[5.4] Gene-level matrix files:"
ls -lh "$MATRIX_DIR"/

echo "  - Trinity.gene.counts.matrix      (gene-level raw counts, for DESeq2)"
echo "  - Trinity.gene.TPM.notMatrix      (gene-level TPM, for visualization)"
echo "  - Trinity.gene.TMM.TPM.matrix     (gene-level TMM-normalized TPM, for WGCNA)"

echo "========================================"
echo "[Step 5] Complete: $(date)"
echo "========================================"
conda deactivate || true
