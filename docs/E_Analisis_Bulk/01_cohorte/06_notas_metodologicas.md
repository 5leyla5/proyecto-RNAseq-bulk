# TCGA-BRCA: definicion preliminar de cohorte

Fuente: GDC Data Portal, Repository (consulta/descarga 2026-11-16).
Proyecto TCGA-BRCA; Gene Expression Quantification; RNA-Seq; STAR - Counts;
Transcriptome Profiling; acceso Open; Tissue Type Tumor; Tumor Descriptor Primary.
Se comprobo adicionalmente codigo TCGA 01 en cada Sample ID.

## Unidad de analisis y elegibilidad
Una unidad clinica = un Case ID; una muestra de expresion por paciente.
En pacientes con varios archivos se priorizo OCT, luego Unknown, luego FFPE;
los empates se resolvieron por Sample ID y File ID. Es una regla preliminar
que debera revisarse tras examinar la calidad y el efecto de lote.

El piloto provisional de 50 pacientes excluye FFPE y utiliza set.seed(20261116).
La seleccion debera revisarse tras obtener datos clinicos (estadio y faltantes).
Las 1111 matrices de expresion NO han sido descargadas por este script.
Los metadatos originales permanecen en data/raw/TCGA_BRCA/ (excluido de Git).

## Comprobaciones reproducibles
Se verifico igualdad por File ID y nombre entre Sample Sheet y Manifest;
tambien se verificaron proyecto, tipo de dato, tipo y descriptor tumoral,
codigo de muestra 01 y coincidencia Case ID / Sample ID.

## Pendiente
Descargar clinical.tsv/clinical.json del GDC, evaluar estadio AJCC y faltantes,
confirmar seleccion piloto, descargar y verificar MD5 de expresion,
comprobar correspondencia de genes con GSE176078 y hacer deconvolucion.
