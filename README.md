# Proyecto RNA-seq bulk | Bioinformática II

**Deconvolución celular de RNA-seq bulk en cáncer de mama con referencia scRNA-seq de GSE176078 y aplicación exploratoria a TCGA-BRCA.**  
**Estado del proyecto:** diseño/factibilidad y preparación de datos documentados; deconvolución **no ejecutada**.  
**Repositorio:** https://github.com/5leyla5/proyecto-RNAseq-bulk  
**Corte de metadatos GDC:** 16-11-2026. **Entrega parcial 1:** 19-11-2026.

## 1. Pregunta, hipótesis y límites

**Pregunta:** ¿qué error presenta una referencia de tipos celulares construida con scRNA-seq al estimar proporciones conocidas en mezclas pseudobulk de donantes independientes y qué tan factible es aplicarla a tumores primarios de TCGA-BRCA?

**Hipótesis verificable:** la deconvolución basada en una referencia de donantes de entrenamiento logrará menor error absoluto medio (MAE) en pseudobulk de donantes reservados que la predicción basal de composición promedio del conjunto de entrenamiento. La aplicación a TCGA-BRCA será **exploratoria**, sin proporciones reales de referencia ni inferencia causal sobre estadio.

## 2. Fuentes y acceso

### scRNA-seq, GSE176078 / Broad SCP1039

- Artículo: Wu, S. Z. et al. (2021), *Nature Genetics*. DOI: https://doi.org/10.1038/s41588-021-00911-1.
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
- Broad Single Cell Portal: https://singlecell.broadinstitute.org/single_cell/study/SCP1039
- Archivos locales originales en `data/raw/GSE176078/`: `Allcells_raw_count_out.zip` y `Whole_miniatlas_meta.csv` (colocados manualmente; no los descarga `scripts/01_descargar_datos.R`).
- Archivo ZIP: matriz de conteos UMI formato 10X dentro de `BrCa_Atlas_Count_out/` (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`).
- Cohorte verificada: **26 donantes**, **100.064 células**, **29.733 genes**; subtipos **11 ER+, 5 HER2+, 10 TNBC**, nueve tipos celulares principales.
- Matriz de origen y comprobantes MD5: consultar `docs/01_registro_datos.csv`. El objeto Seurat/almacén BPCells en `data/processed/GSE176078/` no se publica.

### RNA-seq bulk, TCGA-BRCA / GDC

- Portal de datos: https://portal.gdc.cancer.gov/projects/TCGA-BRCA
- Documentación de cuantificación STAR Counts: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
- Filtros de la exportación (16-11-2026): proyecto `TCGA-BRCA`, `Gene Expression Quantification`, `Transcriptome Profiling`, `RNA-Seq`, `STAR - Counts`, acceso `open`, `Tissue Type = tumor`, `Tumor Descriptor = primary`. Código de muestra TCGA `01` confirmado además en el Sample Sheet.
- Archivos originales colocados en `data/raw/TCGA_BRCA/`: `gdc_sample_sheet.tsv` (1111 filas), `gdc_manifest.txt` (1111), `cohort_TCGA_BRCA_LSQ.2026-11-16.tsv` (1095 UUID de caso); bajo `clinical/` los cinco TSV clínicos de **TCGA-BRCA** (`clinical.tsv`, `exposure.tsv`, `family_history.tsv`, `follow_up.tsv`, `pathology_detail.tsv`); bajo `biospecimen/` los TSV de muestra, porción, alícuota, etc. **No mezclar con la descarga previa de TCGA-CESC.**
- Se conservaron los comprimidos clínicos originales localmente. Solo los archivos de metadatos y resultados **agregados** derivados van a GitHub; datos crudos y por paciente quedan excluidos por `.gitignore`.

## 3. Inventario auditado y decisión de piloto

| Indicador | Resultado confirmado |
|---|---:|
| Archivos STAR Counts en el catálogo | **1.111** |
| Pacientes independientes (Case ID) | **1.095** |
| Identificadores de muestra únicos (Sample ID) | **1.106** |
| Pacientes con múltiples archivos | **11** |
| Pacientes con una muestra no FFPE seleccionable | **1.093** |
| Pacientes con estadio patológico AJCC utilizable | **1.069 / 1.095** |
| Distribución AJCC en cohorte completa | I: 181; II: 619; III: 249; IV: 20; sin dato: 26 |
| Pacientes en piloto inicial | **50** (1 archivo por paciente) |
| Estadios en el piloto | **I: 6; II: 33; III: 10; IV: 1** |
| Preservación de archivos piloto | **OCT: 34; Unknown: 16; FFPE: 0** |
| Tamaño del piloto por manifest | **212.050.127 bytes (~212,05 MB)** |

**Regla del E01:** una muestra por paciente; ordenar por preservación `OCT`, `Unknown`, `FFPE`, luego por `Sample ID` y `File ID`. Excluir FFPE **solo del piloto**; entre los restantes, seleccionar aleatoriamente 50 con `set.seed(20261116)`. **E03 no selecciona en función del estadio ni reemplaza pacientes**: comprueba que los 50 originales tienen estadio utilizable y publica su manifiesto final.

La muestra de 50 sirve para *factibilidad técnica* y para una exploración descriptiva. No está diseñada para una comparación con potencia estadística por estadio: **solo un paciente está en estadio IV**. La preservación `Unknown` no es un método confirmado, por lo que constituye un posible confusor técnico.

## 4. Estructura del proyecto

```text
data/
  raw/GSE176078/                    # insumos originales, excluidos de Git
  raw/TCGA_BRCA/                    # Sample Sheet, Manifest, clínica, biospecimen
  raw/TCGA_BRCA/expresion_piloto/   # matrices descargadas con gdc-client (futuro)
  processed/GSE176078/             # Seurat y BPCells, excluidos de Git
  processed/TCGA_BRCA/             # tablas por paciente, excluidas de Git
scripts/
  01_descargar_datos.R              # VERIFICA insumos preexistentes; NO descarga
  A_Analisis_SingleCell/             # scripts A/01–A/04, QC realizado
  B_Referencia_Celular/             # pendiente
  C_Pseudobulk_Deconvolucion/        # pendiente
  D_Evaluacion_Modelo/               # pendiente
  E_Analisis_Bulk/
    01_preparar_cohorte_TCGA_BRCA.R
    02_validar_clinica_TCGA_BRCA.R
    03_cerrar_piloto_TCGA_BRCA.R     # ejecutar, auditar y versionar
    04_verificar_descarga_TCGA_BRCA.R # ejecutar tras descargar
  F_Discusion_Critica/              # pendiente
  00_generar_informe_entrega1.py   # genera Word y Markdown desde CSV
docs/E_Analisis_Bulk/
  01_cohorte/                        # resultados E01, ya generados
  02_clinica/                        # resultados E02, ya generados
  03_piloto/                         # se genera con E03
  04_descarga/                       # se genera con E04
docs/Entrega_1/                    # informe Word, Markdown y trazabilidad
requirements-informe.txt             # python-docx 1.2.0, separado de R/renv
figures/                             # figuras del proyecto
renv.lock                            # lista de versiones, restauración no validada
README.md
```

## 5. Reproducción en Positron (Windows)

La **consola R** ejecuta `source(...)`; la **terminal PowerShell** ejecuta `git` y `gdc-client`. Abrir la raíz del proyecto (donde están `README.md` y `scripts/`). No se requieren paquetes adicionales para E01–E04: usan R base y `tools`.

### A. Preparar insumos

1. Clonar `https://github.com/5leyla5/proyecto-RNAseq-bulk` y abrir el directorio en Positron.
2. Colocar los datos de GSE176078 en `data/raw/GSE176078/` y los tres archivos de metadata de GDC bajo `data/raw/TCGA_BRCA/`, con los nombres citados. Descomprimir **solo** los TSV clínicos de TCGA-BRCA en `data/raw/TCGA_BRCA/clinical/`. Guardar originales fuera de Git.
3. Ejecutar las etapas A existentes según sus scripts y salidas; consultar las notas y checksums de `docs/`. **No asumir que la reconstrucción limpia de BPCells en otro equipo ya ha sido validada.**

### B. Catalogar y consolidar la cohorte (ya ejecutados en el proyecto de origen)

```r
source('scripts/E_Analisis_Bulk/01_preparar_cohorte_TCGA_BRCA.R')
source('scripts/E_Analisis_Bulk/02_validar_clinica_TCGA_BRCA.R')
```

E01 genera las tablas de `docs/E_Analisis_Bulk/01_cohorte/`, incluido un manifest provisional de 50 archivos. E02 genera `docs/E_Analisis_Bulk/02_clinica/` y `data/processed/TCGA_BRCA/clinica_por_paciente.csv` **(solo local)**. El mensaje final original de E01 sobre ausencia de clínica quedó desactualizado después de que se descargó: no implica que el proceso haya fallado.

### C. Cerrar y auditar el piloto (preparado, pendiente de ejecutar localmente)

```r
source('scripts/E_Analisis_Bulk/03_cerrar_piloto_TCGA_BRCA.R')
```

Esperado: 50 pacientes únicos, 50 estadios disponibles; estadios 6/33/10/1; 34 OCT y 16 Unknown. Se publican `docs/E_Analisis_Bulk/03_piloto/01_resumen_piloto.csv`, `02_estadios_piloto.csv`, `03_manifest_piloto_50_FINAL.txt`, `04_criterios_piloto.md`, `05_md5_insumos.csv` y se conserva el cruce por paciente solo en `data/processed/TCGA_BRCA/`.

### D. Descargar *solo* los 50 archivos de expresión (pendiente)

Instalar **GDC Data Transfer Tool** para Windows desde https://gdc.cancer.gov/access-data/gdc-data-transfer-tool y añadir `gdc-client.exe` a `PATH`, o invocarlo por su ruta completa. En **PowerShell**, desde la raíz del repositorio:

```powershell
# Comprobar que el ejecutable responde:
gdc-client --help

# Descargar ~212 MB (más espacio temporal) y permitir reanudar:
gdc-client download -m "docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt" -d "data/raw/TCGA_BRCA/expresion_piloto"
```

Si `gdc-client` no se reconoce, por ejemplo `& "C:\\ruta\\gdc-client.exe" download -m ... -d ...`. **No utilizar el manifest de 1.111 archivos** para el piloto. No se requiere token para los archivos identificados como acceso abierto. El cliente organiza la descarga en subdirectorios por UUID y soporta reanudación. Documentación oficial: https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/.

### E. Comprobar integridad (solo DESPUÉS de descargar)

```r
source('scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R')
```

E04 compara existencia, tamaño y **MD5** de los 50 archivos contra el manifest. Solo se marca descarga **verificada** si las 50 verificaciones son correctas. Salidas: `docs/E_Analisis_Bulk/04_descarga/` (publicables, sin matrices grandes).

### F. Etapas de análisis pendientes

1. Inspeccionar columnas `gene_id`, `unstranded` y anotaciones de STAR; utilizar **conteos crudos**, no TPM ni FPKM, para procedimientos que requieran conteos. El pipeline GDC cuantifica con STAR sobre GRCh38 y proporciona varias escalas.
2. Conciliar identificadores de genes del bulk y referencia scRNA (Ensembl con/sin versión frente a símbolos, duplicados y genes compartidos). Fijar una única regla de correspondencia antes de generar estimaciones.
3. Separar donantes GSE176078 antes de construir marcadores/perfiles; validar pseudobulk de donantes reservados con MAE/RMSE y comparación basal. No tratar varias mezclas del mismo donante como muestras biológicas independientes.
4. Solo tras validación interna, aplicar referencia al piloto TCGA-BRCA y analizar con prudencia la relación descriptiva de proporciones estimadas y estadio, controlando preservación y factores disponibles; ausencia de verdad celular en bulk real.

## 6. Entorno y reproducibilidad

- Equipo de trabajo: **Windows + Positron**, **R 4.5.3**, **Seurat 5.6.0**, **SeuratObject 5.4.0**, **BPCells 0.3.1**, `renv` **1.3.1** según `renv.lock`/sesión documentada. Revisar lockfile para versiones completas; repositorio base A–F y dependencias del código.
- `renv.lock` existe, pero **no está demostrado que `renv::restore()` funcione desde un clon limpio**, y faltan archivos de activación de renv como `.Rprofile` y `renv/activate.R` (no confundir snapshot con entorno plenamente portable).
- Recomendación para siguiente revisión: instalar/activar renv en una **copia o clon de prueba**, restaurar y documentar problemas antes de modificar el ambiente de trabajo actual.
- La carpeta `data/processed/GSE176078/bpcells_counts/` contiene datos auxiliares usados por Seurat; los RDS pueden depender de rutas o del almacén BPCells externo, por lo que hay que verificar portabilidad con reconstrucción desde el ZIP de origen.
- Para reproducibilidad, conservar la exportación original, su fecha, la semilla, criterios, nombres exactos y hashes. No publicar datos crudos ni rutas de OneDrive personales.

## 7. Resultados publicados frente a trabajos pendientes

**Ejecutado y validado:** acceso y QC scRNA (A/01–A/04), inventario TCGA-BRCA (E01), validación clínica (E02). **Verificado con los archivos originales para preparar E03:** el piloto original contiene 50 casos con estadio y manifiesto de ~212 MB. **Código preparado pero por ejecutar en Positron:** E03 y E04 (este último tras la descarga). **No ejecutados:** descarga de matrices bulk, normalización/UMAP final, B–D, deconvolución real, evaluación cuantitativa del modelo y asociación clínica.

## 8. Regeneración del informe académico y publicación en GitHub

El informe de **tres páginas** se entrega en `docs/Entrega_1/Informe_Entrega_1.docx` y
`docs/Entrega_1/Informe_Entrega_1.md`. **Ambos se construyen desde el mismo script**, que
lee las tablas agregadas `docs/E_Analisis_Bulk/01_cohorte/`, `02_clinica/` y `03_piloto/`;
no incorpora números copiados a mano ni lee matrices de expresión.
La fuente ejecutable del documento es `scripts/00_generar_informe_entrega1.py` y la
dependencia fijada está en `requirements-informe.txt`.

**Orden exacto para reproducir el informe:** ejecutar E01, E02 y E03 en la consola R;
después ejecutar en la **terminal PowerShell** de Positron:

```powershell
py -3 -m pip install -r requirements-informe.txt
py -3 scripts/00_generar_informe_entrega1.py
```

El generador comprueba la coherencia de denominadores, los 50 archivos únicos y los
**212.050.127 bytes** del manifiesto del corte del 16-11-2026. Si los datos de origen
cambian, se debe revisar y versionar de forma explícita el nuevo corte; no se debe
alterar silenciosamente el informe histórico.

La ejecución del generador requiere **Python 3 y python-docx 1.2.0**. Es una dependencia
**separada** del entorno de R/`renv.lock`. El informe está preparado en formato A4 con
saltos de página definidos; la paginación puede variar si un equipo sustituye fuentes.
La versión entregada se verificó en tres páginas; la reconstrucción local debe inspeccionarse
antes de entregar al profesor.

Para verificar y publicar los archivos nuevos desde la raíz del repositorio:

```powershell
git add scripts/E_Analisis_Bulk/03_cerrar_piloto_TCGA_BRCA.R
git add scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R
git add scripts/00_generar_informe_entrega1.py requirements-informe.txt
git add docs/E_Analisis_Bulk/03_piloto docs/Entrega_1 README.md
git diff --cached --check
git status
git commit -m "Entrega 1: piloto TCGA-BRCA e informe reproducible"
git push origin main
```

**Auditoría de origen y resultados:** consultar `docs/Entrega_1/00_trazabilidad_resultados.md`
y `docs/E_Analisis_Bulk/03_piloto/05_md5_insumos.csv`. Los resultados de E01/E02
quedan publicados; E03 debe ejecutarse y confirmar sus salidas antes del `commit`.
No declarar el paso E04 ni la descarga como completados hasta obtener 50 MD5 correctos.

## 9. Referencias fundamentales

1. Wu SZ et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics* 53:1334–1347. https://doi.org/10.1038/s41588-021-00911-1
2. Wang X et al. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. *Nature Communications* 10:380. https://doi.org/10.1038/s41467-018-08023-x
3. Tran KA et al. (2023). Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures. *Nature Communications* 14:5758. https://doi.org/10.1038/s41467-023-41385-5
4. Luecken MD & Theis FJ (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology* 15:e8746. https://doi.org/10.15252/msb.20188746
5. NCI Genomic Data Commons. Guías de Repository, mRNA Analysis Pipeline y Data Transfer Tool. https://docs.gdc.cancer.gov/
