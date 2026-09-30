#!/bin/bash
#=========================================#
# Step 4: TransDecoder ORF prediction     #
#=========================================#
# Predicts ORFs from Trinity transcripts.
# Outputs: .pep (protein), .cds (CDS), .gff3 (annotation).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"

cd "$WORKDIR"

activate_env "$CONDA_RNASEQ" || exit 1

echo "========================================"
echo "[Step 4] TransDecoder ORF Prediction"
echo "Started: $(date)"
echo "========================================"

echo "[4.0] Waiting for Trinity.fasta..."
while [ ! -f "$TRINITY_FA" ]; do
    sleep 60
done
sleep 10
echo "[4.0] Trinity.fasta found: $(date)"

mkdir -p "$TRANSDECODER_DIR"
cd "$TRANSDECODER_DIR"

# ---- 4.1: LongOrfs ----
echo "[4.1] TransDecoder.LongOrfs (min ${MIN_ORF_LENGTH} aa)..."
echo "  Started: $(date)"
TransDecoder.LongOrfs -t ../"$TRINITY_FA" -m "$MIN_ORF_LENGTH"
echo "  Done: $(date)"

# ---- 4.2: Predict ----
echo "[4.2] TransDecoder.Predict..."
echo "  Started: $(date)"
TransDecoder.Predict -t ../"$TRINITY_FA"
echo "  Done: $(date)"

# ---- Summary ----
echo "[4.3] Results:"
echo "  PEP count: $(grep -c '>' Trinity.fasta.pep 2>/dev/null || echo 0)"
echo "  CDS count: $(grep -c '>' Trinity.fasta.cds 2>/dev/null || echo 0)"
echo "  GFF lines: $(wc -l < Trinity.fasta.gff3 2>/dev/null || echo 0)"
ls -lh Trinity.fasta.pep Trinity.fasta.cds Trinity.fasta.gff3 2>/dev/null || true

echo "========================================"
echo "[Step 4] Complete: $(date)"
echo "========================================"
conda deactivate || true
