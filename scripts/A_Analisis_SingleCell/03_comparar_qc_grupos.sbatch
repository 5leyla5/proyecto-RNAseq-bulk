#!/bin/bash
#SBATCH --job-name=bioII_A03
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/A_Analisis_SingleCell/03_qc_por_grupos/slurm_%j.out
#SBATCH --error=logs/A_Analisis_SingleCell/03_qc_por_grupos/slurm_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
unset R_HOME

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

R_BIN="/mnt/beegfs/apps/r-4.5.1/bin/Rscript"

echo "=== OBJETIVO A / ETAPA 03 ==="
echo "Nodo: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Inicio: $(date)"

"$R_BIN" --vanilla -e '
options(bitmapType = "cairo")
source("scripts/A_Analisis_SingleCell/03_comparar_qc_grupos.R")
'

echo "Finalizacion: $(date)"
