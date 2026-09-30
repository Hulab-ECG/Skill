#!/bin/bash
#=========================================#
# Step 3: Reference genome + HISAT2 index #
#=========================================#
# Reference genome FASTA is required.
#
# Annotation auto-detect:
#   - GTF and/or GFF present: decompress, if needed convert GFF -> GTF with
#     gffread -T, extract reference transcripts/proteins, and build HISAT2
#     index with the annotation.
#   - Neither GTF nor GFF present: build a plain HISAT2 index.  StringTie
#     assembly + merge + TransDecoder will be handled in Step 5.

set -euo pipefail
export ZSH_VERSION=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"
cd "$WORKDIR"

activate_env "$CONDA_RNASEQ"

echo "========================================"
echo "[Step 3] Reference Genome + HISAT2 Index"
echo "Started: $(date)"
echo "========================================"

mkdir -p "$REF_DIR"
mkdir -p "$ANNOTATION_DIR"

# ---- 3.1: Reference genome FASTA (must be supplied locally) ----
echo "[3.1] Checking reference genome: $REF_ASSEMBLY"
if [ ! -f "$REF_GENOME" ] && [ -f "${REF_DIR}/${REF_GENOME_NAME}.gz" ]; then
    echo "  Decompressing existing ${REF_GENOME_NAME}.gz ..."
    gunzip -f "${REF_DIR}/${REF_GENOME_NAME}.gz"
fi
if [ ! -s "$REF_GENOME" ]; then
    echo "  ERROR: reference genome not found: $REF_GENOME" >&2
    echo "  Place your genomic FASTA in $REF_DIR as $REF_GENOME (or ${REF_GENOME_NAME}.gz)" >&2
    exit 934
fi
echo "  Genome exists: $(grep -c '>' "$REF_GENOME") contigs"
ls -lh "$REF_GENOME"

# ---- 3.2: Reference GTF annotation (optional) ----
echo "[3.2] Checking reference GTF annotation: $REF_GTF"
if [ ! -f "$REF_GTF" ] && [ -f "${REF_DIR}/${REF_GTF_GZ}" ]; then
    echo "  Decompressing existing ${REF_GTF_GZ} ..."
    gunzip -f "${REF_DIR}/${REF_GTF_GZ}"
fi
if [ -s "$REF_GTF" ]; then
    echo "  GTF: $REF_GTF ($(grep -c $'\ttranscript\t' "$REF_GTF" 2>/dev/null || echo 0) transcripts)"
    ls -lh "$REF_GTF"
else
    echo "  No GTF found locally (${REF_GTF} or ${REF_GTF_GZ})"
fi

# ---- 3.3: Reference GFF annotation (optional) ----
echo "[3.3] Checking reference GFF annotation: $REF_GFF"
if [ ! -f "$REF_GFF" ] && [ -f "${REF_DIR}/${REF_GFF_GZ}" ]; then
    echo "  Decompressing existing ${REF_GFF_GZ} ..."
    gunzip -f "${REF_DIR}/${REF_GFF_GZ}"
fi
if [ -s "$REF_GFF" ]; then
    echo "  GFF: $REF_GFF"
else
    echo "  No GFF found locally (${REF_GFF} or ${REF_GFF_GZ})"
fi

# Make sure StringTie GTF exists; default expectation is that GTF is
# already supplied, otherwise convert from GFF (StringTie -G requires GTF).
if [ -s "$REF_GFF" ] || [ -s "$REF_GTF" ]; then
    echo "[3.3b] Ensuring StringTie GTF exists..."
    if [ ! -s "$REF_GTF" ]; then
        echo "  Converting GFF -> GTF with gffread -T ..."
        gffread "$REF_GFF" -T -o "$REF_GTF"
    fi
    echo "  GTF ready: $REF_GTF ($(grep -c $'\ttranscript\t' "$REF_GTF" 2>/dev/null || echo 0) transcripts)"

    # Extract reference transcript/protein sequences with gffread.
    echo "[3.4] Extracting reference transcripts/proteins with gffread..."
    REFPEP_GFFREAD="${ANNOTATION_DIR}/ref_transcripts.pep"
    REFCDS_GFFREAD="${ANNOTATION_DIR}/ref_transcripts.fasta"
    if [ ! -s "$REFCDS_GFFREAD" ]; then
        gffread "$REF_GFF" -g "$REF_GENOME" -w "$REFCDS_GFFREAD"
    fi
    if [ ! -s "$REFPEP_GFFREAD" ]; then
        gffread "$REF_GFF" -g "$REF_GENOME" -y "$REFPEP_GFFREAD"
    fi
    echo "  Transcripts: $(grep -c '>' "$REFCDS_GFFREAD" || true)"
    echo "  Proteins:    $(grep -c '>' "$REFPEP_GFFREAD" || true)"

    # Build transcript/protein -> gene ID maps used by downstream analyses.
    echo "[3.4b] Building reference annotation ID maps..."
    python3 "$SCRIPT_DIR/reference_prot_gene_map.py" \
        --gff "$REF_GFF" \
        --gtf "$REF_GTF" \
        --out "$PROT_TO_GENE_MAP" \
        --trans_to_gene "$TRANSCRIPT_TO_GENE_MAP"
else
    echo "[3.3b] No GTF or GFF found."
    echo "  Step 5 will use StringTie assembly + merge + TransDecoder to build the transcript annotation."
fi

# ---- 3.5: Build HISAT2 index ----
echo "[3.5] Building HISAT2 index..."
if [ -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "  Index already exists"
else
    hisat2-build -p "$THREADS_HISAT2" "$REF_GENOME" "$HISAT2_INDEX"
fi
ls -lh "${REF_DIR}"/hisat2_index.*

echo "========================================"
echo "[Step 3] Complete: $(date)"
echo "========================================"
conda deactivate
