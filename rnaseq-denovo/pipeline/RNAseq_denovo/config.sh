#=========================================#
# RNAseq_denovo Pipeline                  #
# De Novo RNA-seq quant:
#   fastp QC -> Trinity assembly       ->
#   TransDecoder ORF prediction        ->
#   Salmon quant -> count/TPM matrices
#
# Steps 1-2 are shared with the reference pipeline:
#   - 01_download_rawdata.sh
#   - 02_fastp_qc.sh
# Steps 3-5 follow the de novo portable workflow:
#   - 03_trinity_assembly.sh
#   - 04_transdecoder.sh
#   - 05_salmon_quant.sh
#=========================================#

PROJECT_NAME="RNAseq_denovo"
BIOPROJECT=""   # example: PRJNA642665; set to your BioProject or leave empty

# ---- Portability ----
# CONDA_BASE: leave empty to auto-detect from $CONDA_PREFIX or common paths.
# WORKDIR: leave empty to use the directory that contains this pipeline.
CONDA_BASE=""
WORKDIR=""

# ---- Conda environments ----
CONDA_RNASEQ="rnaseq"   # Trinity, TransDecoder, Salmon, HMMER, seqkit, sra-tools, fastp

# ---- Samples ----
# Put one sample SRA accession per line, e.g.
#   SRR12640679
#   SRR12640678
# Leave empty or set to your own runs.
CONTROL_SAMPLES=""

# ---- Threads & Memory ----
THREADS_DOWNLOAD=10
THREADS_QC=10
THREADS_SORT=10
THREADS_TRINITY=40
THREADS_SALMON=10
MAX_MEMORY="100G"

# ---- Trinity parameters ----
TRINITY_NORMALIZE="--normalize_reads"
TRINITY_EXTRA="--full_cleanup"

# ---- fastp parameters ----
FASTP_PARAMS="--cut_tail -W 4 -M 20 -n 10 -c -l 50"

# ---- TransDecoder ----
MIN_ORF_LENGTH=100

# ---- Salmon ----
SALMON_KMER=31
SALMON_TYPE="quasi"

# ---- Paths (relative to WORKDIR) ----
RAWDATA_DIR="00.rawdata"
CLEANDATA_DIR="01.cleandata"
SORTED_DIR="01.cleandata.sorted"
TRINITY_DIR="02.trinity/trinity_out_dir"
TRANSDECODER_DIR="03.transdecoder"
SALMON_DIR="04.salmon"
MATRIX_DIR="04.salmon/matrix"
ANNOTATION_DIR="05.annotation"

# ---- Key output files ----
TRINITY_FA="${TRINITY_DIR}/Trinity.fasta"
GENE_TRANS_MAP="${TRINITY_DIR}/Trinity.fasta.gene_trans_map"
PEP_FILE="${TRANSDECODER_DIR}/Trinity.fasta.pep"
CDS_FILE="${TRANSDECODER_DIR}/Trinity.fasta.cds"
COUNT_MATRIX="${MATRIX_DIR}/Trinity.gene.counts.matrix"
TPM_MATRIX="${MATRIX_DIR}/Trinity.gene.TMM.TPM.matrix"

echo "[config] ${PROJECT_NAME} (de novo, Trinity + TransDecoder -> Salmon)"
echo "[config] Workdir: ${WORKDIR:-<auto>}"
echo "[config] BioProject: ${BIOPROJECT}"
echo "[config] Threads (Trinity): ${THREADS_TRINITY}"
echo "[config] Max memory: ${MAX_MEMORY}"

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
    if command -v TransDecoder.LongOrfs >/dev/null 2>&1; then
        cand="$(dirname "$(command -v TransDecoder.LongOrfs)")/../opt/transdecoder/TransDecoder"
        if [ -x "$cand" ]; then
            printf '%s\n' "$cand"
            return 0
        fi
    fi
    printf '\n'
}
