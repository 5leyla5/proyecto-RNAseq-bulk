#!/bin/bash
#SBATCH --job-name=bioII_A04
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/A_Analisis_SingleCell/04_filtrado_qc/slurm_%j.out
#SBATCH --error=logs/A_Analisis_SingleCell/04_filtrado_qc/slurm_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
unset R_HOME

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

R_BIN="/mnt/beegfs/apps/r-4.5.1/bin/Rscript"

echo "=== OBJETIVO A / ETAPA 04 ==="
echo "Nodo: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Inicio: $(date)"

"$R_BIN" --vanilla \
  scripts/A_Analisis_SingleCell/04_documentar_decision_qc.R

echo "Finalizacion: $(date)"
EOFcd ~/proyectos/proyecto-RNAseq-bulk

mkdir -p logs/A_Analisis_SingleCell/04_filtrado_qc

cat > scripts/A_Analisis_SingleCell/04_documentar_decision_qc.sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=bioII_A04
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/A_Analisis_SingleCell/04_filtrado_qc/slurm_%j.out
#SBATCH --error=logs/A_Analisis_SingleCell/04_filtrado_qc/slurm_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
unset R_HOME

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

R_BIN="/mnt/beegfs/apps/r-4.5.1/bin/Rscript"

echo "=== OBJETIVO A / ETAPA 04 ==="
echo "Nodo: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Inicio: $(date)"

"$R_BIN" --vanilla \
  scripts/A_Analisis_SingleCell/04_documentar_decision_qc.R

echo "Finalizacion: $(date)"
