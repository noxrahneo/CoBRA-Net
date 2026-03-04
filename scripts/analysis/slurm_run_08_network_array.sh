#!/usr/bin/env bash
set -euo pipefail

# Slurm array launcher for network building + visualization (08c + 08d)
# Usage examples:
#   sbatch scripts/analysis/slurm_run_08_network_array.sh
#   sbatch --array=0-5 scripts/analysis/slurm_run_08_network_array.sh
#
# Optional env vars:
#   CONDA_ENV=breast_cancer_scrnaseq
#   NETWORK_TYPE=signed
#   MAX_EDGES=2500
#   MIN_WEIGHT=0.10
#   POWERS=1,2,3,4,5,6,7,8,9,10,12,14,16,18,20
#   INTERACTIVE_HTML=1

#SBATCH --job-name=net08
#SBATCH --output=logs/slurm/%x_%A_%a.out
#SBATCH --error=logs/slurm/%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --array=0-5

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

mkdir -p logs/slurm

CONDA_ENV="${CONDA_ENV:-breast_cancer_scrnaseq}"
NETWORK_TYPE="${NETWORK_TYPE:-signed}"
MAX_EDGES="${MAX_EDGES:-2500}"
MIN_WEIGHT="${MIN_WEIGHT:-0.10}"
POWERS="${POWERS:-1,2,3,4,5,6,7,8,9,10,12,14,16,18,20}"
INTERACTIVE_HTML="${INTERACTIVE_HTML:-1}"

CONDITIONS=(
  "ER_tumor"
  "HER2_tumor"
  "Normal"
  "Normal_BRCA1_-_pre-neoplastic"
  "Triple_negative_BRCA1_tumor"
  "Triple_negative_tumor"
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$IDX" -lt 0 || "$IDX" -ge "${#CONDITIONS[@]}" ]]; then
  echo "Invalid array index: $IDX"
  exit 1
fi

CONDITION="${CONDITIONS[$IDX]}"
echo "[$(date)] Condition=$CONDITION"

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate "$CONDA_ENV"
else
  echo "conda not found on PATH"
  exit 1
fi

python3 scripts/analysis/08c_network_power_tom_prep.py \
  --condition "$CONDITION" \
  --input-dir results/stages/07_network/correlation/pearson \
  --output-dir results/stages/07_network/network_prep \
  --network-type "$NETWORK_TYPE" \
  --powers "$POWERS"

VIS_ARGS=(
  --condition "$CONDITION"
  --input-dir results/stages/07_network/network_prep
  --output-dir results/stages/07_network/network_viz
  --network-type "$NETWORK_TYPE"
  --max-edges "$MAX_EDGES"
  --min-weight "$MIN_WEIGHT"
)

if [[ "$INTERACTIVE_HTML" == "1" ]]; then
  VIS_ARGS+=(--interactive-html)
fi

python3 scripts/analysis/08d_networkx_visualization.py "${VIS_ARGS[@]}"

echo "[$(date)] Finished Condition=$CONDITION"
