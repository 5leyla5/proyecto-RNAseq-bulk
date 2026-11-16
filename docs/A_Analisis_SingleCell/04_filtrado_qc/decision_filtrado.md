# Objetivo A / 04 - Decision documentada de control de calidad

## Decision
Se conservaron las celulas y los genes del archivo scRNA-seq GSE176078 recibido
del Broad Single Cell Portal, sin aplicar filtros adicionales en esta etapa.
Los conteos crudos y los objetos previos NO fueron modificados.
Esta es una decision fundamentada de QC, no la afirmacion de que todas
las celulas son perfectas ni de que se descarta revisar posibles dobletes.

## Punto de partida
- Celulas conservadas: **100.064** (eliminadas adicionales: **0**).
- Genes conservados: **29.733** (eliminados adicionales: **0**).
- Donantes: **26**.
- Tipos celulares principales: **9**.
- Capa utilizada: conteos UMI crudos del ensayo RNA, respaldados por BPCells.

## Fundamento del control de calidad celular
Wu et al. (2021) procesaron los datos originales con Cell Ranger y EmptyDrops,
e indicaron >200 genes detectados, >250 UMI y <20 % mitocondrial.
El objeto recibido ya satisface estos tres criterios en todas sus celulas:
- nFeature_RNA > 200: TRUE.
- nCount_RNA > 250: TRUE.
- percent.mt < 20: TRUE.
En A/02 se recalculo el porcentaje mitocondrial desde los conteos crudos,
contrastandolo con los metadatos; en A/03 se compararon los indicadores
por tipo celular y por donante.
Las distribuciones son heterogeneas; un unico umbral adicional muy estricto
podria excluir de manera desproporcionada poblaciones epiteliales tumorales
u otras poblaciones de interes para la referencia de deconvolucion.

## Evaluacion exploratoria del QC de genes
Se definio un gen como detectado en una celula si el conteo UMI es >= 1.
Los escenarios siguientes fueron simulados en A/03, pero NO aplicados:

| Deteccion minima | Genes conservados | Genes excluidos hipoteticos | Porcentaje excluido |
|---|---:|---:|---:|
| >=1 celulas | 29.733 | 0 | 0,00 % |
| >=3 celulas | 27.719 | 2.014 | 6,77 % |
| >=10 celulas | 25.209 | 4.524 | 15,22 % |

Se evaluaron ademas la distribucion de frecuencia de deteccion, la relacion
entre UMI medio crudo y deteccion, y los 20 genes con mayor UMI acumulado.
La baja deteccion global no basta para descartar un gen: podria ser relevante
en una poblacion celular poco abundante.
La seleccion de genes altamente variables (HVG) para PCA sera una etapa
posterior y diferente del filtrado de genes del objeto original.

## Celulas conservadas por tipo celular
| Tipo celular | Celulas conservadas |
|---|---:|
| B-cells | 3.206 |
| CAFs | 6.573 |
| Cancer Epithelial | 24.489 |
| Endothelial | 7.605 |
| Myeloid | 9.675 |
| Normal Epithelial | 4.355 |
| Plasmablasts | 3.524 |
| PVL | 5.423 |
| T-cells | 35.214 |

Los recuentos por donante y por combinacion donante-tipo celular estan
documentados en `celulas_por_donante.csv` y `celulas_por_donante_y_tipo.csv`.
Incluyen las combinaciones sin celulas y no representan proporciones reales
de tejido tumoral, sino distribucion de celulas capturadas y anotadas.

## Limitaciones y decisiones futuras
- El atlas ya tiene un preprocesamiento previo; no estamos evaluando gotas vacias originales.
- El QC descriptivo no descarta por si solo dobletes, contaminacion ambiental o efectos de lote.
- No se seleccionaran marcadores de referencia con datos de donantes de evaluacion.
- En el objetivo B se separaran donantes antes de seleccionar los marcadores de la referencia.
- Para pseudobulk se conservaran y emplearan conteos crudos UMI, no la capa normalizada.
- Si aparecieran problemas posteriores, cualquier nuevo filtro se documentara
  con nuevos parametros, recuentos por donante y tipo celular y una nueva version del analisis.

## Reproducibilidad
- Ejecutar desde la raiz del proyecto:
  `source("scripts/A_Analisis_SingleCell/04_documentar_decision_qc.R")`.
- Objeto de entrada: `data/processed/GSE176078/seurat_qc_inicial.rds`.
- Resumen de genes de entrada: `docs/A_Analisis_SingleCell/03_qc_por_grupos/qc_genes_impacto_umbrales_pre.csv`.
- Este script solo lee resultados de A/02 y A/03 y escribe tablas y este documento.
- Fuente de los archivos originales y sumas MD5: `docs/01_registro_datos.csv`.
- Versiones de software: `logs/A_Analisis_SingleCell/04_filtrado_qc/session_info.txt`.

## Referencias
1. Wu SZ et al. (2021). A single-cell and spatially resolved atlas of human breast cancers.
   Nature Genetics. https://doi.org/10.1038/s41588-021-00911-1
2. Single-cell best practices: Quality Control.
   https://www.sc-best-practices.org/preprocessing-visualization/quality-control/

**Decision registrada: no aplicar filtros adicionales a celulas ni genes.**
