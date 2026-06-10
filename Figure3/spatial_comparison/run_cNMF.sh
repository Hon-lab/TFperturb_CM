#!/usr/bin/env bash
#SBATCH -J cnmf_run
#SBATCH -p 256GBv2
#SBATCH -t 99-1:00:00
#SBATCH -o logs/log_%x.out
#SBATCH -e logs/log_%x.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

#To use cNMF
source activate geneprogram
target_config_path="${SCRIPT_DIR}/config_cnmf_spatial.json"

echo "config=${target_config_path}"

python "${SCRIPT_DIR}/cnmf_run.py" "${target_config_path}"
