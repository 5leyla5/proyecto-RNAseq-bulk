# E/05 - Lectura y QC inicial de expresion bulk TCGA-BRCA

**Estado:** resultados derivados de las 50 matrices STAR Counts descargadas y verificadas por E/04.
Este script NO produce deconvolucion ni valida biologicamente proporciones celulares.

## Entradas y procedencia

- Cohorte: manifiesto FINAL E/03, corte GDC 2026-11-16, 50 archivos tumor primario.
- Integridad: se exigen 50 MD5 y tamanos correctos; E/05 vuelve a calcular MD5.
- Columnas GDC: `gene_id`, `gene_name`, `gene_type`, `unstranded`.
- `unstranded` corresponde a conteos crudos. No se usan TPM ni FPKM.
- Se excluyen de la matriz genica las cuatro filas tecnicas de STAR: `N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`.
- No hay filtro de genes ni muestras en E/05; las metricas son descriptivas.

## Resultados observados

- Archivos bulk leidos: **50**.
- Genes cuantificados por muestra: **60660**.
- Genes con conteos >0 en al menos una de las 50 muestras: **52937**.
- Mediana de conteos genicos totales por muestra: **55487415**.
- Mediana de genes detectados por muestra: **32180**.
- Referencia GSE176078: **Leida desde ZIP GSE176078**.
- Coincidencias exactas por simbolo, no ambiguas: **20288**.

## Interpretacion y limitaciones

La matriz local conserva los ID Ensembl con VERSION; no se sustituyen por simbolos.
GSE176078 en 10x presenta nombres de genes; la coincidencia por simbolo es
una evaluacion preliminar, no una armonizacion validada ni un mapeo uno-a-uno
de IDs. Genes con nombres duplicados se excluyen de la comparacion por simbolo,
pero NO de la matriz de conteos. Antes de B/C/E final se necesita resolver
versiones, sinonimos, anotaciones y genes comunes de forma explicita.
El total de conteos genicos NO representa el total de lecturas secuenciadas.
No se elimina ninguna muestra por QC en esta etapa; revisar distribuciones
y predefinir reglas antes de filtrar. No hay estimacion de celulas en bulk.

## Ubicacion de resultados

- Tablas agregadas para GitHub: `docs/E_Analisis_Bulk/05_expresion/`.
- Graficos para GitHub: `figures/E_Analisis_Bulk/05_expresion/`.
- Conteos, anotaciones y metricas individuales solo locales: `data/processed/TCGA_BRCA/05_expresion/`.
- Entorno: `05_sessionInfo.txt`; codigo fuente: `scripts/E_Analisis_Bulk/05_leer_QC_expresion_TCGA_BRCA.R`.

## Documentacion fuente

- GDC mRNA Pipeline: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
