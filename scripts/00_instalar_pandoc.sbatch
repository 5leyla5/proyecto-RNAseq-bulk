#!/bin/bash
#SBATCH --job-name=bioII_PANDOC
#SBATCH --partition=ngen-ko
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/Entrega_1/pandoc_%j.out
#SBATCH --error=logs/Entrega_1/pandoc_%j.err

set -euo pipefail

echo "=== INSTALACION PANDOC ==="
echo "Nodo: $(hostname)"
echo "Inicio: $(date)"

BASE="$HOME/.local/opt"
DESCARGAS="$HOME/.local/downloads"
VERSION="3.14"
ARCHIVO="pandoc-3.14-linux-amd64.tar.gz"
BIN="$BASE/pandoc-$VERSION/bin"

SHA256="c5f4720cdefb3f7bdc5d07202ef579feef0d1370773c873804dd8fdaa8c4f1b8"

mkdir -p "$BASE" "$DESCARGAS"

if [ ! -x "$BIN/pandoc" ]; then
    cd "$DESCARGAS"

    if [ ! -f "$ARCHIVO" ]; then
        curl -fL --retry 3 --retry-delay 5 \
          --connect-timeout 20 --max-time 900 \
          -o "$ARCHIVO" \
          "https://github.com/jgm/pandoc/releases/download/$VERSION/$ARCHIVO"
    fi

    echo "$SHA256  $ARCHIVO" | sha256sum -c -

    tar -xzf "$ARCHIVO" -C "$BASE"
fi

test -x "$BIN/pandoc"

export PATH="$BIN:$PATH"
export RSTUDIO_PANDOC="$BIN"
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.5"

echo "=== VERSION PANDOC ==="
pandoc --version | head -n 2

echo "=== VERIFICACION EN R ==="

/mnt/beegfs/apps/r-4.5.1/bin/Rscript --vanilla -e '
stopifnot(rmarkdown::pandoc_available())
cat("Pandoc disponible:", rmarkdown::pandoc_available(), "\n")
cat("Version:", as.character(rmarkdown::pandoc_version()), "\n")
'

echo "PANDOC INSTALADO CORRECTAMENTE"
echo "Finalizacion: $(date)"
