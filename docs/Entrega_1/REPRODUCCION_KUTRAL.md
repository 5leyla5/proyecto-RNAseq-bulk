# Reproducción de la Entrega 1 en Kutral

## 1. Entorno de ejecución

Clúster Kutral, Universidad de Chile.
Partición Slurm: ngen-ko.

El nodo de acceso host04 se utiliza únicamente para
editar scripts, preparar archivos, enviar trabajos a
Slurm y consultar resultados.

Los análisis R se ejecutan en nodos de cómputo.

## 2. Versiones utilizadas

- R: 4.5.1
- Seurat: 6.0.0
- SeuratObject: 5.4.0
- BPCells: 0.3.1
- Matrix: 1.7.3
- Pandoc para el informe: 3.14

Ejecutable R:
/mnt/beegfs/apps/r-4.5.1/bin/Rscript

Biblioteca personal:
~/R/x86_64-pc-linux-gnu-library/4.5

El archivo renv.lock registra el entorno de Windows,
con R 4.5.3 y Seurat 5.6.0. No representa exactamente
el entorno utilizado en Kutral.

La restauración completa desde un clon limpio todavía
no ha sido validada.

## 3. Preparación de registros Slurm

Desde la raíz del repositorio:

    bash scripts/00_preparar_logs_slurm.sh

Este paso crea los directorios que Slurm necesita
antes de iniciar los trabajos.

Los registros permanecen en logs/, fuera de GitHub.

## 4. Insumos originales

GSE176078:
- Allcells_raw_count_out.zip
- Whole_miniatlas_meta.csv

Ambos deben estar en data/raw/GSE176078/.

Sus tamaños y MD5 están documentados en:
docs/01_registro_datos.csv

La reconstrucción de los objetos Seurat y BPCells
requiere acceso a estos archivos originales.

## 5. Ejecución del objetivo A

Enviar los trabajos en el siguiente orden,
comprobando que cada uno finalice correctamente:

    sbatch scripts/A_Analisis_SingleCell/01_crear_objeto_seurat.sbatch
    sbatch scripts/A_Analisis_SingleCell/02_qc_inicial.sbatch
    sbatch scripts/A_Analisis_SingleCell/03_comparar_qc_grupos.sbatch
    sbatch scripts/A_Analisis_SingleCell/04_documentar_decision_qc.sbatch

No enviar el siguiente trabajo hasta verificar
que el anterior terminó con COMPLETED y ExitCode 0:0.

No se realizan normalización, PCA, UMAP ni
deconvolución en esta etapa de la entrega.

## 6. Ejecuciones verificadas

- A/01: trabajo 1103125, COMPLETED.
- A/02: trabajo 1103127, COMPLETED.
- A/03: trabajo 1103128, COMPLETED.
- A/04: trabajo 1103481, COMPLETED.
- Informe: trabajo 1103450, COMPLETED.

El trabajo A/04 anterior, 1103130, falló debido
a instrucciones accidentales dentro del archivo
Slurm. El script fue corregido y ejecutado
nuevamente sin errores.

## 7. Regeneración del informe

Primero preparar los directorios de registros.

Luego ejecutar:

    sbatch scripts/00_render_informe_entrega1.sbatch

Salida:
docs/Entrega_1/Informe_Entrega_1.html

El informe se genera a partir de CSV y manifiestos
versionados. No requiere descargar las matrices
originales ni recalcular los análisis.

## 8. TCGA-BRCA

La descarga de 1093 archivos y su validación MD5
se completaron en Windows.

Evidencia:
docs/E_Analisis_Bulk/06_cohorte_completa/integridad_1093.csv

Script de verificación:
scripts/E_Analisis_Bulk/08_validar_integridad_1093.ps1

Los 1093 archivos originales todavía no se han
transferido a Kutral.

El QC de expresión de la cohorte completa
permanece pendiente. Los resultados E/05
corresponden únicamente al piloto de 50 muestras.

## 9. Limitaciones de reproducibilidad

- Los datos originales no están versionados en Git.
- Los objetos BPCells requieren sus archivos en disco.
- La instalación de paquetes en un entorno limpio
  no se ha comprobado completamente.
- El entorno Windows difiere del entorno Kutral.
- La ejecución de todos los pasos desde un clon
  nuevo todavía debe someterse a prueba.

Este documento registra ejecuciones verificadas,
pero no declara reproducibilidad integral ya probada.
