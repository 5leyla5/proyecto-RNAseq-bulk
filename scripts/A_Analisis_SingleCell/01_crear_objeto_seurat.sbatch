#!/bin/bash
#SBATCH --job-name=bioII_A01
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/A_Analisis_SingleCell/01_objeto_seurat/slurm_%j.out
#SBATCH --error=logs/A_Analisis_SingleCell/01_objeto_seurat/slurm_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
unset R_HOME

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1

R_BIN="/mnt/beegfs/apps/r-4.5.1/bin/Rscript"

echo "=== OBJETIVO A / ETAPA 01 ==="
echo "Nodo: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Inicio: $(date)"

"$R_BIN" --vanilla \
  scripts/A_Analisis_SingleCell/01_crear_objeto_seurat.R

echo "Finalizacion: $(date)"
