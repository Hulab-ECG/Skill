#!/bin/bash
#=========================================#
# RNAseq_denovo Pipeline                  #
#   raw FASTQ check/download -> fastp QC  #
#   -> Trinity assembly -> TransDecoder   #
#   -> Salmon quant -> count/TPM matrices #
#=========================================#
# Usage:
#   bash run_pipeline.sh           # Run all steps
#   bash run_pipeline.sh [1|2|3|4|5] # Run a single step

set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
source "${PIPELINE_DIR}/config.sh"
cd "$WORKDIR"

run_step() {
    local script=$1
    echo ""
    echo "############################################################"
    echo "## Running $(basename "$script")"
    echo "############################################################"
    bash "$script"
}

if [ $# -eq 0 ]; then
    run_step "$SCRIPTS_DIR/01_download_rawdata.sh"
    run_step "$SCRIPTS_DIR/02_fastp_qc.sh"
    run_step "$SCRIPTS_DIR/03_trinity_assembly.sh"
    run_step "$SCRIPTS_DIR/04_transdecoder.sh"
    run_step "$SCRIPTS_DIR/05_salmon_quant.sh"

    echo ""
    echo "============================================================"
    echo "  DE NOVO QUANTIFICATION PIPELINE COMPLETE!"
    echo "  $(date)"
    echo "============================================================"
    echo "  Transcripts:   $TRINITY_FA"
    echo "  Proteins:      $PEP_FILE"
    echo "  CDS:           $CDS_FILE"
    echo "  Count matrix:  $COUNT_MATRIX"
    echo "  TMM TPM matrix: $TPM_MATRIX"
    echo "============================================================"
else
    case $1 in
        1) run_step "$SCRIPTS_DIR/01_download_rawdata.sh" ;;
        2) run_step "$SCRIPTS_DIR/02_fastp_qc.sh" ;;
        3) run_step "$SCRIPTS_DIR/03_trinity_assembly.sh" ;;
        4) run_step "$SCRIPTS_DIR/04_transdecoder.sh" ;;
        5) run_step "$SCRIPTS_DIR/05_salmon_quant.sh" ;;
        *)
            echo "Unknown step: $1 (available: 1, 2, 3, 4, 5)"
            exit 955
            ;;
    esac
fi
