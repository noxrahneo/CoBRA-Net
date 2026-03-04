#!/bin/bash
#SBATCH --job-name=net08_all
#SBATCH --output=/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net/logs/net08_%j.out
#SBATCH --error=/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net/logs/net08_%j.err
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --nodelist=antony

set -euo pipefail

PROJECT_ROOT="/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net"
SCRIPT_DIR="$PROJECT_ROOT/scripts/analysis"

CORR_IN="$PROJECT_ROOT/results/stages/07_network/correlation/pearson"
PREP_OUT="$PROJECT_ROOT/results/stages/07_network/network_prep"
VIZ_OUT="$PROJECT_ROOT/results/stages/07_network/network_viz"

echo "[$(date)] SLURM job ${SLURM_JOB_ID}: running 08c + 08d (all conditions)"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate breast_cancer_scrnaseq

cd "$SCRIPT_DIR"

python3 08c_network_power_tom_prep.py \
  --condition all \
  --input-dir "$CORR_IN" \
  --output-dir "$PREP_OUT" \
  --network-type signed

python3 08d_networkx_visualization.py \
  --condition all \
  --input-dir "$PREP_OUT" \
  --output-dir "$VIZ_OUT" \
  --network-type signed \
  --max-edges 1400 \
  --min-weight 0.14 \
  --label-top-hubs 20 \
  --hub-neighbors 30 \
  --interactive-html \
  --interactive-max-edges 900 \
  --interactive-min-weight 0.20

echo "[$(date)] Finished 08c + 08d"
