#=========================================#
# RNAseq_reference Pipeline               #
# Reference-based RNA-seq quant:
#   fastp QC -> reference/+ HISAT2 index ->
#   HISAT2 align -> StringTie quant ->
#   count + TPM matrices
#
# Annotation mode is auto-detected:
#   - GTF/GFF present  -> StringTie -eB -G ref.gtf
#   - neither present  -> StringTie assembly + merge + TransDecoder
#=========================================#

PROJECT_NAME="RNAseq_reference"
BIOPROJECT=""   # example: PRJNA642665; set to your BioProject or leave empty

# ---- Portability ----
# CONDA_BASE: leave empty to auto-detect from $CONDA_PREFIX or common paths.
# WORKDIR: leave empty to use the directory that contains this pipeline.
CONDA_BASE=""
WORKDIR=""

# ---- Conda environments ----
CONDA_RNASEQ="rnaseq"   # stringtie, hisat2, samtools, gffread, seqkit, sra-tools, fastp, transdecoder

# Optional: path to TransDecoder wrapper if it is NOT inside $CONDA_RNASEQ/bin.
# Leave empty for auto-detection (checks inside the RNA-seq env first).
TRANS_DECODER_WRAP=""

# ---- Samples ----
# Put one sample SRA accession per line, e.g.
#   SRR12640679
#   SRR12640678
# Leave empty or set to your own runs.
CONTROL_SAMPLES=""

# ---- Reference genome + annotations (supplied locally) ----
# Example values:
#   REF_ASSEMBLY="GCA_014636745.1"
#   REF_GENOME_NAME="GCA_014636745.1_ASM1463674v1_genomic.fna"
#   REF_GTF_GZ="GCA_014636745.1_ASM1463674v1_genomic.gtf.gz"
#   REF_GFF_GZ="GCA_014636745.1_ASM1463674v1_genomic.gff.gz"
REF_ASSEMBLY=""
REF_GENOME_NAME=""
REF_GTF_GZ=""
REF_GFF_GZ=""

# ---- Threads ----
THREADS_DOWNLOAD=10
THREADS_QC=10
THREADS_HISAT2=40
THREADS_STRINGTIE=40

# ---- fastp ----
FASTP_PARAMS="--cut_tail -W 4 -M 20 -n 10 -c -l 50"

# ---- Paths ----
RAWDATA_DIR="00.rawdata"
CLEANDATA_DIR="01.cleandata"
REF_DIR="02.ref_genome"
ALIGN_DIR="03.hisat2"
STRINGTIE_DIR="04.stringtie"
COUNT_MATRIX="${STRINGTIE_DIR}/gene_count_matrix.csv"
TPM_MATRIX="${STRINGTIE_DIR}/gene_TPM_matrix.tsv"
ANNOTATION_DIR="05.annotation"

REF_GENOME="${REF_DIR}/${REF_GENOME_NAME}"
REF_GTF="${REF_DIR}/${REF_GTF_GZ%.gz}"
REF_GFF="${REF_DIR}/${REF_GFF_GZ%.gz}"
HISAT2_INDEX="${REF_DIR}/hisat2_index"

# Annotation outputs
#   GTF/GFF mode:          gffread transcripts/proteins
#   no-annotation mode:    gffread transcripts + TransDecoder proteins/CDS
MERGED_GTF="${STRINGTIE_DIR}/merged.gtf"
TRANSCRIPTS_FASTA="${ANNOTATION_DIR}/transcripts.fasta"
REF_PEP="${ANNOTATION_DIR}/transcripts.fasta.transdecoder.pep"
REF_CDS="${ANNOTATION_DIR}/transcripts.fasta.transdecoder.cds"
PROT_TO_GENE_MAP="${REF_DIR}/protein_to_gene_id.txt"
TRANSCRIPT_TO_GENE_MAP="${REF_DIR}/transcript_to_gene_id.txt"

echo "[config] ${PROJECT_NAME} (reference-based quant, auto annotation mode)"
echo "[config] Reference: ${REF_ASSEMBLY}"
echo "[config] Main steps: fastp QC -> HISAT2 -> StringTie quant -> count/TPM matrix"
echo "[config] If no GTF/GFF is placed in ${REF_DIR}, Step 5 will run StringTie assembly + merge + TransDecoder."

# ================================================================
#  AUTO-DETECTION — do not edit below this line
# ================================================================
_PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-${0}}")" && pwd)"

if [ -z "${WORKDIR:-}" ]; then
    WORKDIR="$_PIPELINE_ROOT"
fi
mkdir -p "$WORKDIR"

if [ -z "${CONDA_BASE:-}" ]; then
    if [ -n "${CONDA_PREFIX:-}" ]; then
        CONDA_BASE="$(dirname "$CONDA_PREFIX")"
    else
        for cand in "$HOME/miniconda3" "$HOME/anaconda3" "/opt/miniconda3" "/opt/anaconda3"; do
            if [ -d "$cand" ]; then
                CONDA_BASE="$cand"
                break
            fi
        done
    fi
fi
if [ -z "${CONDA_BASE:-}" ]; then
    echo "[config] ERROR: conda base not found. Set CONDA_BASE in config.sh" >&2
fi

if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    CONDA_INIT_SCRIPT="${CONDA_BASE}/etc/profile.d/conda.sh"
else
    CONDA_INIT_SCRIPT=""
fi

activate_env() {
    local env_name="$1"
    if [ -z "$CONDA_INIT_SCRIPT" ] || [ ! -f "$CONDA_INIT_SCRIPT" ]; then
        echo "[activate_env] ERROR: conda.sh not found ($CONDA_INIT_SCRIPT)" >&2
        return 998
    fi
    # shellcheck disable=SC1090
    source "$CONDA_INIT_SCRIPT" || return 099
    conda activate "$env_name" || return 009
    export PATH="${CONDA_BASE}/envs/${env_name}/bin:$PATH"
}

control_sample_list() {
    printf '%s\n' "$CONTROL_SAMPLES" \
        | sed -e 's/[[:space:]]//g' -e '/^$/d'
}

# Locate the TransDecoder wrapper that takes `-t transcripts.fasta`.
find_transdecoder() {
    if [ -n "${TRANS_DECODER_WRAP:-}" ] && [ -x "${TRANS_DECODER_WRAP:-}" ]; then
        printf '%s\n' "$TRANS_DECODER_WRAP"
        return 0
    fi
    local cand longorfs_base cand2
    for cand in \
        "${CONDA_BASE}/envs/${CONDA_RNASEQ}/opt/transdecoder/TransDecoder" \
        "${CONDA_BASE}/envs/${CONDA_RNASEQ}/bin/TransDecoder" \
        ; do
        if [ -n "$cand" ] && [ -x "$cand" ]; then
            printf '%s\n' "$cand"
            return 0
        fi
    done
    longorfs_base="$(command -v TransDecoder.LongOrfs 2>/dev/null || true)"
    if [ -n "$longorfs_base" ]; then
        cand2="$(dirname "$longorfs_base")/../opt/transdecoder/TransDecoder"
        if [ -x "$cand2" ]; then
            printf '%s\n' "$cand2"
            return 0
        fi
    fi
    printf '\n'
}
