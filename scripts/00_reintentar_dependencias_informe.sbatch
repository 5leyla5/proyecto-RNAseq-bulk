#!/bin/bash
#SBATCH --job-name=bioII_DEP2
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/Entrega_1/dependencias_%j.out
#SBATCH --error=logs/Entrega_1/dependencias_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
unset R_HOME

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1

/mnt/beegfs/apps/r-4.5.1/bin/Rscript --vanilla - <<'RSCRIPT'

options(timeout = 600)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

paquetes <- c(
  "bslib",
  "rmarkdown",
  "knitr",
  "base64enc"
)

cat("=== PAQUETES ANTES DE INSTALAR ===\n")

for (p in paquetes) {
  cat(p, requireNamespace(p, quietly = TRUE), "\n")
}

# Instalar solo los paquetes faltantes.
for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) {

    cat("\nInstalando:", p, "\n")

    install.packages(
      p,
      repos = "https://cloud.r-project.org",
      lib = Sys.getenv("R_LIBS_USER"),
      Ncpus = 1
    )

    if (!requireNamespace(p, quietly = TRUE)) {
      stop("No se pudo instalar: ", p)
    }
  }
}

cat("\n=== VERSIONES FINALES ===\n")

for (p in paquetes) {
  cat(p, as.character(packageVersion(p)), "\n")
}

cat("\n=== PANDOC ===\n")

cat("Ruta en PATH:", Sys.which("pandoc"), "\n")
cat("Disponible para R Markdown:",
    rmarkdown::pandoc_available(), "\n")

if (rmarkdown::pandoc_available()) {
  cat("Version:",
      as.character(rmarkdown::pandoc_version()), "\n")
}

cat("\nDEPENDENCIAS R COMPLETADAS\n")

RSCRIPT
