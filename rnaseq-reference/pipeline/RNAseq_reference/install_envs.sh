#!/bin/bash
set -euo pipefail

# One-shot deployment helper:
#   1. create/update conda environments from envs/*.yml
#   2. verify required command-line tools are available
# After running this, edit config.sh (CONDA_BASE/WORKDIR) and run run_pipeline.sh.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v conda >/dev/null 2>&1; then
    echo "[install_envs] ERROR: conda not found in PATH." >&2
    echo "[install_envs] Install Miniconda first: https://docs.conda.io/en/latest/miniconda.html" >&2
    exit 1
fi

CONDA_BASE="${CONDA_BASE:-}"
if [ -z "$CONDA_BASE" ]; then
    CONDA_BASE="$(dirname "$(dirname "$(command -v conda)")")"
fi

# shellcheck disable=SC1090
source "${CONDA_BASE}/etc/profile.d/conda.sh"

create_or_update_env() {
    local name="$1"
    local yaml="$2"
    if conda env list | grep -qE "^${name}[[:space:]]"; then
        echo "[install_envs] Updating conda env '${name}' from ${yaml} ..."
        conda env update -f "$yaml" --prune
    else
        echo "[install_envs] Creating conda env '${name}' from ${yaml} ..."
        conda env create -f "$yaml"
    fi
}

create_or_update_env "rnaseq" "envs/env_rnaseq.yml"

echo ""
echo "[install_envs] Verifying command availability ..."

conda activate rnaseq
for cmd_ in stringtie hisat2 samtools gffread seqkit fastp TransDecoder.LongOrfs TransDecoder.Predict; do
    if command -v "$cmd_" >/dev/null 2>&1; then
        echo "  OK  ${cmd_}"
    else
        echo "  MISSING  ${cmd_}"
    fi
done
conda deactivate

echo ""
echo "[install_envs] Done."
echo "[install_envs] Next:"
echo "  1. Edit config.sh and set CONDA_BASE / WORKDIR if needed."
echo "  2. Run: bash run_pipeline.sh"
