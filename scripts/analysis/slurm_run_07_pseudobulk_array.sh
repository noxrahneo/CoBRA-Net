#!/bin/bash
#SBATCH --job-name=pb07_all
#SBATCH --output=/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net/logs/pb07_%j.out
#SBATCH --error=/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net/logs/pb07_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --nodelist=antony

set -euo pipefail

PROJECT_ROOT="/triumvirate/home/alexarol/breast_cancer_analysis/CoBRA-Net"
SCRIPT_DIR="$PROJECT_ROOT/scripts/analysis"

echo "[$(date)] SLURM job ${SLURM_JOB_ID}: running all conditions sequentially on antony"

# Activate conda env
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate breast_cancer_scrnaseq

cd "$SCRIPT_DIR"

python3 07_pseudobulk.py \
  --condition all \
  --threshold-mode auto \
  --input-dir "$PROJECT_ROOT/results/stages/04_annotation_rdata" \
  --output-dir "$PROJECT_ROOT/results/stages/07_network/pseudobulk"

echo "[$(date)] Finished all conditions"
