#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p \
  logs/A_Analisis_SingleCell/01_objeto_seurat \
  logs/A_Analisis_SingleCell/02_qc_inicial \
  logs/A_Analisis_SingleCell/03_qc_por_grupos \
  logs/A_Analisis_SingleCell/04_filtrado_qc \
  logs/Entrega_1

echo "Directorios de Slurm preparados correctamente."
