#!/bin/bash
#=========================================#
# Step 5: StringTie quantification        #
#=========================================#
# Annotation mode is auto-detected:
#  A) Reference GTF/GFF exists -> StringTie -eB -G ref.gtf (no assembly)
#  B) No GTF/GFF exists        -> StringTie assembly + merge -> gffread
#                                -> TransDecoder (PEP/CDS) -> StringTie -eB
#
# Both modes finish with:
#   prepDE.py -> gene_count_matrix.csv
#   build_tpm_matrix.py -> gene_TPM_matrix.tsv

set -euo pipefail
export ZSH_VERSION=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$(dirname "$SCRIPT_DIR")/config.sh"
cd "$WORKDIR"

activate_env "$CONDA_RNASEQ"

echo "========================================"
echo "[Step 5] StringTie Quantification"
echo "Started: $(date)"
echo "========================================"

if [ -s "$REF_GTF" ] || [ -s "$REF_GFF" ]; then
    MODE="reference"
    if [ ! -s "$REF_GTF" ] && [ -s "$REF_GFF" ]; then
        echo "[Step 5] Converting GFF -> GTF with gffread -T ..."
        mkdir -p "$REF_DIR"
        gffread "$REF_GFF" -T -o "$REF_GTF"
    fi
    echo "[Step 5] Using reference GTF: $REF_GTF"
else
    MODE="assembly"
    echo "[Step 5] No GTF/GFF detected; using StringTie assembly + merge + TransDecoder."
fi

STRINGTIE=$(which stringtie)
PREPDE="${SCRIPT_DIR}/prepDE.py"

mkdir -p "$STRINGTIE_DIR"
mkdir -p "$ANNOTATION_DIR"

# One control sample ID -> sorted BAM path (empty if missing)
bam_for_sample() {
    local s="$1"
    local b="${ALIGN_DIR}/${s}.sorted.bam"
    if [ -f "$b" ]; then printf '%s\n' "$b"; fi
}

# Wait for sorted BAM files of the control samples to finish alignment.
echo "[5.0] Waiting for sorted BAM files (control samples only)..."
while [ -z "$(for s in $(control_sample_list); do bam_for_sample "$s"; done)" ]; do
    sleep 60
done
sleep 5
echo "[5.0] BAM files ready: $(for s in $(control_sample_list); do bam_for_sample "$s"; done | wc -l)"

if [ "$MODE" = "reference" ]; then

    # ---- 5.1a: StringTie -eB against reference GTF ----
    echo "[5.1] StringTie -eB quantification with $REF_GTF ..."
    for SampleID in $(control_sample_list); do
        bam="$(bam_for_sample "$SampleID")"
        if [ -z "$bam" ]; then
            echo "  WARNING: BAM missing for $SampleID, skipping"
            continue
        fi
        SAMPLE_DIR="${STRINGTIE_DIR}/${SampleID}_quant"
        mkdir -p "$SAMPLE_DIR"
        out_gtf="${SAMPLE_DIR}/${SampleID}.transcripts.gtf"
        if [ -f "${SAMPLE_DIR}/${SampleID}.gene_abund.tab" ] && [ -f "$out_gtf" ]; then
            echo "  Skipping $SampleID (already quantified)"
            continue
        fi
        echo "  >>> $SampleID ($(date))"
        $STRINGTIE "$bam" \
            -p "$THREADS_STRINGTIE" \
            -G "$REF_GTF" \
            -eB \
            -o "$out_gtf" \
            -A "${SAMPLE_DIR}/${SampleID}.gene_abund.tab"
    done
    QUANT_GTF_SUFFIX=".transcripts.gtf"

else

    # ---- 5.1b: StringTie assembly per sample ----
    echo "[5.1] StringTie assembly per sample..."
    for SampleID in $(control_sample_list); do
        bam="$(bam_for_sample "$SampleID")"
        if [ -z "$bam" ]; then
            echo "  WARNING: BAM missing for $SampleID, skipping assembly"
            continue
        fi
        out_gtf="${STRINGTIE_DIR}/${SampleID}.gtf"
        if [ -f "$out_gtf" ]; then
            echo "  Skipping $SampleID (already assembled)"
            continue
        fi
        echo "  >>> $SampleID ($(date))"
        $STRINGTIE "$bam" -p "$THREADS_STRINGTIE" -o "$out_gtf"
    done

    # ---- 5.2b: Merge transcripts ----
    echo "[5.2] Merging transcripts (stringtie --merge)..."
    if [ ! -f "$MERGED_GTF" ]; then
        : > "${STRINGTIE_DIR}/gtf_list.txt"
        for SampleID in $(control_sample_list); do
            g="${STRINGTIE_DIR}/${SampleID}.gtf"
            [ -f "$g" ] && echo "$g" >> "${STRINGTIE_DIR}/gtf_list.txt"
        done
        if [ -s "${STRINGTIE_DIR}/gtf_list.txt" ]; then
            $STRINGTIE --merge -p "$THREADS_STRINGTIE" \
                -o "$MERGED_GTF" \
                $(cat "${STRINGTIE_DIR}/gtf_list.txt")
        else
            echo "  ERROR: no per-sample StringTie GTF files found." >&2
            exit 1
        fi
    fi
    echo "  Merged transcripts: $(grep -c $'\ttranscript\t' "$MERGED_GTF")"
    echo "  Merged genes:      $(grep -c $'\tgene\t' "$MERGED_GTF")"

    # ---- 5.3b: Extract transcript sequences ----
    echo "[5.3] Extracting transcript sequences (gffread)..."
    GFFREAD=$(which gffread || true)
    if [ -z "$GFFREAD" ]; then
        echo "  ERROR: gffread not found." >&2
        exit 1
    fi
    if [ ! -f "$TRANSCRIPTS_FASTA" ]; then
        $GFFREAD "$MERGED_GTF" -g "$REF_GENOME" -w "$TRANSCRIPTS_FASTA"
    fi
    echo "  Transcripts: $(grep -c '>' "$TRANSCRIPTS_FASTA")"

    # ---- 5.4b: TransDecoder (protein prediction) ----
    TDWRAP="$(find_transdecoder)"
    if [ -z "$TDWRAP" ]; then
        echo "[5.4] ERROR: TransDecoder wrapper not found. Install TransDecoder or set TRANS_DECODER_WRAP in config.sh." >&2
        exit 1
    fi
    echo "[5.4] TransDecoder.LongOrfs + Predict..."
    if [ ! -f "$REF_PEP" ]; then
        cd "$ANNOTATION_DIR"
        $TDWRAP -t transcripts.fasta
        cd "$WORKDIR"
    fi
    echo "  Predicted proteins: $(grep -c '>' "$REF_PEP" 2>/dev/null || echo 0)"

    # ---- 5.5b: StringTie -eB using merged.gtf ----
    echo "[5.5] StringTie -eB quantification with $MERGED_GTF ..."
    for SampleID in $(control_sample_list); do
        bam="$(bam_for_sample "$SampleID")"
        if [ -z "$bam" ]; then
            echo "  WARNING: BAM missing for $SampleID, skipping quantification"
            continue
        fi
        SAMPLE_DIR="${STRINGTIE_DIR}/${SampleID}_quant"
        mkdir -p "$SAMPLE_DIR"
        out_gtf="${SAMPLE_DIR}/${SampleID}.gtf"
        if [ -f "${SAMPLE_DIR}/${SampleID}.gene_abund.tab" ] && [ -f "$out_gtf" ]; then
            echo "  Skipping $SampleID (already quantified)"
            continue
        fi
        echo "  >>> $SampleID ($(date))"
        $STRINGTIE "$bam" \
            -e -B \
            -G "$MERGED_GTF" \
            -p "$THREADS_STRINGTIE" \
            -o "$out_gtf" \
            -A "${SAMPLE_DIR}/${SampleID}.gene_abund.tab"
    done
    QUANT_GTF_SUFFIX=".gtf"

    # ---- 5.6b: Build transcript -> gene mapping ----
    echo "[5.6] Building transcript -> gene mapping from merged.gtf..."
    python3 "$SCRIPT_DIR/build_prot_gene_map.py" "$MERGED_GTF" "$PROT_TO_GENE_MAP"

fi

# ---- 5.2: prepDE.py -> gene-level count matrix ----
echo "[5.2] Running prepDE.py for gene-level counts..."
SAMPLE_LIST="${STRINGTIE_DIR}/sample_list.txt"
: > "$SAMPLE_LIST"
for sname in $(control_sample_list); do
    dir="${STRINGTIE_DIR}/${sname}_quant/"
    [ -d "$dir" ] || continue
    echo -e "${sname}\t${dir}${sname}${QUANT_GTF_SUFFIX}" >> "$SAMPLE_LIST"
done
cat "$SAMPLE_LIST"

python3 "$PREPDE" \
    -i "$SAMPLE_LIST" \
    -g "$COUNT_MATRIX" \
    -t "$STRINGTIE_DIR/transcript_count_matrix.csv"
echo "  Count matrix: $(wc -l < "$COUNT_MATRIX") genes"

# ---- 5.3/5.7: build gene-level TPM matrix from -A abundance tables ----
echo "[5.3/5.7] Building gene-level TPM matrix (from StringTie -A gene abundance tables)..."
python3 "$SCRIPT_DIR/build_tpm_matrix.py" \
    --abund_dir "$STRINGTIE_DIR" \
    --pattern "gene_abund.tab" \
    --samples "$(control_sample_list | tr '\n' ',')" \
    --output "$TPM_MATRIX"
echo "  TPM matrix: $(wc -l < "$TPM_MATRIX") lines"

echo "========================================"
echo "[Step 5] Complete: $(date)"
echo "  Count matrix: $COUNT_MATRIX"
echo "  TPM matrix:   $TPM_MATRIX"
if [ "$MODE" = "assembly" ]; then
    echo "  Proteins:     $REF_PEP"
    echo "  CDS:          $REF_CDS"
fi
echo "========================================"
conda deactivate
