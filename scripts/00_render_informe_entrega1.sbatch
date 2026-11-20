#!/bin/bash
#SBATCH --job-name=bioII_INF1
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/Entrega_1/slurm_%j.out
#SBATCH --error=logs/Entrega_1/slurm_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

unset R_HOME
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export PATH="$HOME/.local/opt/pandoc-3.14/bin:$PATH"
export RSTUDIO_PANDOC="$HOME/.local/opt/pandoc-3.14/bin"
R_BIN="/mnt/beegfs/apps/r-4.5.1/bin/Rscript"

echo "=== BIOINFORMATICA II / ENTREGA 1 ==="
echo "Nodo: $(hostname)"
echo "Trabajo: $SLURM_JOB_ID"
echo "Inicio: $(date)"

"$R_BIN" --vanilla -e '
options(bitmapType = "cairo")
paquetes <- c("rmarkdown", "knitr", "base64enc")
disponibles <- vapply(
  paquetes, requireNamespace,
  logical(1), quietly = TRUE
)

print(data.frame(
  paquete = paquetes,
  instalado = disponibles
))

if (!all(disponibles)) {
  stop("Faltan dependencias para generar el informe.")
}

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc no disponible en el nodo de computo.")
}

cat("Pandoc:", as.character(rmarkdown::pandoc_version()), "\n")

source(
  "scripts/00_render_informe_entrega1.R",
  encoding = "UTF-8"
)
'

echo "Finalizacion: $(date)"
