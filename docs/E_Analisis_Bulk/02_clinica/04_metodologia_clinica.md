# TCGA-BRCA: disponibilidad de datos clinicos

Fuente: GDC, exportacion de TCGA-BRCA del 2026-11-16.
Entradas locales: data/raw/TCGA_BRCA/clinical/clinical.tsv y
data/raw/TCGA_BRCA/cohort_TCGA_BRCA_LSQ.2026-11-16.tsv.

## Metodo reproducible
- Verificacion estricta de proyecto TCGA-BRCA y concordancia de Case ID.
- Multiples filas de clinical.tsv se consolidan por paciente.
- Estadio AJCC obtenido SOLO de diagnosticos con diagnosis_is_primary_disease=true.
- Stage X, sin dato y ausencia de diagnostico primario: no utilizables.
- Estadios contradictorios del diagnostico primario se marcan como conflicto.
- Se reportan estadios I-IV; la agrupacion conserva un apartado para 0 u otros.
- Los datos de paciente por identificador se dejan en data/processed/ (no GitHub).

## Alcance
Este script valida METADATOS clinicos. No descarga matrices STAR Counts,
no realiza deconvolucion ni afirma asociacion entre estadios y proporciones.
Antes de analizar RNA-seq bulk debe descargarse la subcohorte seleccionada
y comprobarse la correspondencia de genes con la referencia scRNA-seq.
