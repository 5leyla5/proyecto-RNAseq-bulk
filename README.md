# Proyecto RNA-seq bulk | Bioinformática II

**Diseño y factibilidad de un flujo reproducible para estimar composición celular de tumores de mama mediante RNA-seq bulk, usando una referencia scRNA-seq.**

**Entrega Parcial 1:** 19 de noviembre de 2026  
**Corte de metadatos GDC:** 16 de noviembre de 2026  
**Estado de los análisis al 17 de noviembre de 2026; documentación revisada el 18 de noviembre:** objetivo A (acceso y QC inicial) y E/01–E/05 (cohorte, clínica, piloto, descarga verificada y QC descriptivo bulk) ejecutados. **La deconvolución celular no se ha ejecutado.**  
**Responsable:** Leyla Victoria Saavedra Quinteros · Estudiante del Magíster en Informática Médica, Universidad de Chile  
**Repositorio:** https://github.com/5leyla5/proyecto-RNAseq-bulk

> **Alcance de la Entrega 1.** Se documentan selección y acceso a datos, definición de cohortes, unidades de análisis, tabla preliminar de muestras, entorno computacional y factibilidad de lectura de expresión. La separación de donantes para validación, la construcción de la referencia, el pseudobulk, MuSiC/NNLS y sus resultados corresponden a etapas posteriores.

## 1. Pregunta, hipótesis y objetivos

**Pregunta de investigación.** ¿Con qué precisión puede una referencia celular construida a partir de scRNA-seq estimar proporciones conocidas de tipos celulares en mezclas pseudobulk de donantes independientes, y es técnicamente factible aplicar ese procedimiento a tumores primarios de TCGA-BRCA?

**Hipótesis que se evaluará:** un método de deconvolución construido a partir de donantes de entrenamiento obtendrá menor error absoluto medio (MAE) en pseudobulk de donantes reservados que un comparador basal que prediga la composición promedio del entrenamiento. **La hipótesis todavía no ha sido contrastada.**

**Objetivo general.** Desarrollar y evaluar un procedimiento reproducible de estimación de composición celular en RNA-seq bulk de cáncer de mama utilizando scRNA-seq como referencia.

**Objetivos específicos:**

1. Caracterizar y comprobar la factibilidad técnica de GSE176078 (scRNA-seq) y TCGA-BRCA (RNA-seq bulk).
2. Construir una referencia celular con donantes de entrenamiento, sin utilizar información de los donantes reservados para seleccionar genes o ajustar el modelo.
3. Generar pseudobulk con proporciones conocidas a partir de donantes reservados y comparar NNLS y MuSiC, si sus supuestos y dependencias resultan adecuados.
4. Evaluar MAE, RMSE y errores por categoría celular frente al comparador basal.
5. Tras la validación interna, aplicar el procedimiento de forma exploratoria a TCGA-BRCA y discutir sus limitaciones.

## 2. Fuentes, archivos originales y unidades de análisis

### 2.1 GSE176078 / Broad SCP1039 — referencia scRNA-seq

- **Origen:** [GEO GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078), [Broad Single Cell Portal SCP1039](https://singlecell.broadinstitute.org/single_cell/study/SCP1039); atlas descrito por Wu et al. (2021).
- **Insumos descargados manualmente:** `data/raw/GSE176078/Allcells_raw_count_out.zip` y `data/raw/GSE176078/Whole_miniatlas_meta.csv`. Los nombres, tamaños y MD5 de los originales están registrados en [`docs/01_registro_datos.csv`](docs/01_registro_datos.csv).
- El ZIP contiene archivos de la matriz 10X (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`). `scripts/01_descargar_datos.R` **comprueba insumos existentes**, no descarga automáticamente desde GEO o Broad.
- **Inventario inicial:** 26 donantes (11 ER+, 5 HER2+, 10 TNBC), **100.064 células**, **29.733 genes** y **9 categorías celulares principales**.
- **Cobertura desigual por donante:** *Normal Epithelial* se observa en **13/26** donantes, *Plasmablasts* en **14/26** y *Cancer Epithelial* en **20/26**. Se calculó desde [`celulas_por_donante_y_tipo.csv`](docs/A_Analisis_SingleCell/04_filtrado_qc/celulas_por_donante_y_tipo.csv), contando donantes con `n_celulas_conservadas > 0`. Los ceros no prueban ausencia biológica; pueden reflejar captura y anotación. La futura partición B deberá revisar la cobertura del conjunto de entrenamiento.
- **Unidad biológica de independencia:** el *donante*, no cada célula. Las células son observaciones anidadas dentro de sus donantes. La partición entrenamiento/prueba todavía no se ha ejecutado.
- **Estado:** carga, comprobaciones iniciales y QC descriptivo realizados en objetivo A. No interpretar esto como validación final de normalización, anotaciones o modelo.

### 2.2 TCGA-BRCA / Genomic Data Commons — RNA-seq bulk

- **Origen:** [NCI GDC — TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA); [documentación de cuantificación STAR Counts](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).
- **Exportación de catálogo GDC (16-11-2026):** proyecto `TCGA-BRCA`; categoría `Transcriptome Profiling`; tipo `Gene Expression Quantification`; estrategia `RNA-Seq`; flujo `STAR - Counts`; acceso `open`; tejido `Tumor`, descriptor `Primary`; código TCGA de muestra `01` comprobado en el Sample Sheet.
- **Entradas locales:** `data/raw/TCGA_BRCA/gdc_sample_sheet.tsv`, `gdc_manifest.txt` y `cohort_TCGA_BRCA_LSQ.2026-11-16.tsv`; archivos clínicos correctos de **TCGA-BRCA** en `data/raw/TCGA_BRCA/clinical/` y archivos de biospecimen en `data/raw/TCGA_BRCA/biospecimen/`. No mezclar con una exportación previa de `TCGA-CESC`.
- **Unidad clínica:** un `Case ID` = un paciente. **Unidad de expresión:** una muestra/archivo STAR Counts seleccionado por paciente. Un paciente con varios archivos no se contabiliza varias veces como individuo independiente.
- **Inventario:** 1.111 archivos STAR Counts, 1.106 `Sample ID` distintos y **1.095 pacientes** distintos; 11 pacientes tienen múltiples archivos (16 archivos adicionales frente al número de pacientes).
- **Estado:** metadatos y clínica validados; 50 matrices del piloto descargadas e íntegras; lectura y QC descriptivo E/05 ejecutados. Las 1.111 matrices de la cohorte completa **no** se han descargado.

## 3. Cohortes y tabla preliminar para la Entrega 1

| Variable | GSE176078 (scRNA-seq) | TCGA-BRCA (bulk) |
|---|---:|---:|
| Donantes o pacientes identificados | 26 | 1.095 |
| Células individuales | 100.064 | No observadas por bulk |
| Genes de la fuente inicial / anotados | 29.733 | 60.660 |
| Categorías celulares conocidas | 9 principales | No observadas; **por estimar** |
| Archivos de expresión | ZIP con conteos unicelulares | 1.111 STAR Counts catalogados |
| Unidad principal para independencia | Donante | Paciente |
| Datos utilizados hasta esta entrega | Referencia scRNA y QC inicial | Piloto de 50 muestras (1 por paciente) |
| Control inicial | Realizado (objetivo A) | 50/50 MD5 correctos y QC descriptivo E/05 |

### 3.1 Cohorte bulk y disponibilidad clínica

Los scripts [E/01](scripts/E_Analisis_Bulk/01_preparar_cohorte_TCGA_BRCA.R) y [E/02](scripts/E_Analisis_Bulk/02_validar_clinica_TCGA_BRCA.R) generaron los resultados siguientes:

| Indicador | Valor |
|---|---:|
| Pacientes totales con expresión catalogada | **1.095** |
| Pacientes con una muestra no FFPE seleccionable | **1.093** |
| Pacientes con diagnóstico identificado como primario | **1.093** |
| Pacientes con estadio patológico AJCC utilizable | **1.069** |
| Pacientes sin estadio AJCC utilizable | **26** |
| Estadios AJCC I / II / III / IV (cohorte completa) | **181 / 619 / 249 / 20** |
| Pacientes con edad disponible | **1.094** |

**Importante:** el valor 1.093 de elegibilidad no FFPE y el valor 1.093 de diagnóstico primario corresponden a verificaciones **diferentes**; su coincidencia numérica no demuestra que sean exactamente los mismos pacientes. Los datos clínicos originales contienen múltiples filas por paciente y se consolidan antes de producir los resúmenes.

Evidencias: [`01_resumen_TCGA_BRCA.csv`](docs/E_Analisis_Bulk/01_cohorte/01_resumen_TCGA_BRCA.csv), [`01_resumen_clinico.csv`](docs/E_Analisis_Bulk/02_clinica/01_resumen_clinico.csv), [`02_estadios_ajcc.csv`](docs/E_Analisis_Bulk/02_clinica/02_estadios_ajcc.csv).

### 3.2 Selección del piloto de 50 pacientes

1. Seleccionar un archivo tumoral primario por `Case ID`. En caso de varios archivos, priorizar preservación **OCT → Unknown → FFPE**; resolver empates por `Sample ID` y `File ID` en orden determinista.
2. Excluir del **piloto** las muestras FFPE seleccionadas: quedan **1.093** candidatos no FFPE. La categoría `Unknown` indica preservación no informada, no tejido OCT confirmado.
3. Seleccionar aleatoriamente 50 pacientes con **`set.seed(20261116)`** en E/01. El piloto no se estratifica ni se elige con base en estadio.
4. E/03 comprueba la coincidencia exacta de los 50 `File ID` con E/01 y el estadio de los pacientes, **sin sustituir casos**.

**Piloto comprobado:** 50 pacientes / 50 archivos de tumor primario; **34 OCT, 16 Unknown y 0 FFPE**; estadio **I: 6, II: 33, III: 10, IV: 1**. El piloto sirve para **factibilidad técnica**, no tiene potencia para inferencias por estadio. La ampliación futura a hasta 1.093 pacientes no FFPE dependerá del QC, concordancia clínica y disponibilidad de recursos.

Evidencias: [`01_resumen_piloto.csv`](docs/E_Analisis_Bulk/03_piloto/01_resumen_piloto.csv), [`02_estadios_piloto.csv`](docs/E_Analisis_Bulk/03_piloto/02_estadios_piloto.csv), [`03_manifest_piloto_50_FINAL.txt`](docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt), [`04_criterios_piloto.md`](docs/E_Analisis_Bulk/03_piloto/04_criterios_piloto.md).

## 4. Resultados verificados de descarga y QC de expresión (E/04–E/05)

### 4.1 Transferencia de archivos: E/04

Se utilizó **GDC Data Transfer Tool (`gdc-client`) versión 2.3** en PowerShell. Los 50 archivos del manifiesto final suman **212.050.127 bytes** (212,05 MB decimales). La verificación R comparó existencia, tamaño y MD5 contra el manifiesto:

| Métrica E/04 | Resultado |
|---|---:|
| Archivos esperados / presentes | 50 / 50 |
| Archivos con tamaño y MD5 válidos | **50** |
| Ausentes o con error | **0** |
| Estado | **Verificada** |

Evidencias: [`02_resumen_descarga.csv`](docs/E_Analisis_Bulk/04_descarga/02_resumen_descarga.csv), [`01_integridad_archivos.csv`](docs/E_Analisis_Bulk/04_descarga/01_integridad_archivos.csv) y [protocolo reproducible](docs/E_Analisis_Bulk/04_descarga/04_protocolo_descarga.md).

> **Nota sobre archivos históricos:** algunas notas metodológicas generadas en E/01–E/03 conservan expresiones como «pendiente de descargar» porque documentan el estado de cada etapa cuando fueron creadas. El estado **actual**, respaldado por `02_resumen_descarga.csv` y E/05, es **descarga íntegra y lectura de 50 matrices completadas**. No se han alterado retroactivamente esas salidas de scripts.

**Procedencia:** la descarga efectiva se realizó con comandos manuales de PowerShell y `gdc-client.exe`; el script [`04_descargar_bulk_TCGA_BRCA.ps1`](scripts/E_Analisis_Bulk/04_descargar_bulk_TCGA_BRCA.ps1) fue añadido **después** como alternativa automatizada para futuras reproducciones. No se afirma que ese `.ps1` haya sido probado en el equipo de origen.

### 4.2 Lectura y QC inicial de expresión: E/05

El script [`05_leer_QC_expresion_TCGA_BRCA.R`](scripts/E_Analisis_Bulk/05_leer_QC_expresion_TCGA_BRCA.R) volvió a comprobar los MD5 y leyó las 50 matrices con **`gene_id`, `gene_name`, `gene_type` y conteos crudos `unstranded`**; excluyó de la matriz génica las filas técnicas `N_unmapped`, `N_multimapping`, `N_noFeature` y `N_ambiguous`. No utilizó TPM/FPKM, no filtró muestras ni genes y **no ejecutó deconvolución**.

| Indicador E/05 | Resultado observado |
|---|---:|
| Matrices bulk leídas | **50** |
| Genes anotados por matriz | **60.660** |
| Genes `protein_coding` | **19.962** |
| Genes con conteo positivo en ≥1 muestra | **52.937** |
| Genes con conteo positivo en las 50 muestras | **17.680** |
| Modelo génico informado | **GENCODE v36** |
| Mediana de conteos génicos por muestra | **55.487.414,5** |
| Mediana de genes detectados por muestra | **32.180** |
| Símbolos génicos compartidos, no ambiguos | **20.288** |

**Correspondencia génica:** TCGA-BRCA utiliza identificadores Ensembl **con versión** (`gene_id`); el archivo `features.tsv` de GSE176078 aporta símbolos/nombres y no confirma identificadores ENSG equivalentes. Los **20.288 símbolos compartidos** son una comparación descriptiva por coincidencia exacta y no ambigua dentro de cada origen (aprox. **68,2 %** de los 29.733 símbolos de referencia). **No** constituyen todavía una armonización válida para combinar matrices o estimar composición celular. El solapamiento de Ensembl no pudo calcularse con los identificadores disponibles.

**Interpretación del QC:** existe variabilidad en conteos y genes detectados entre muestras, pero el E/05 es **descriptivo**, sin exclusión por umbrales ni prueba de efectos de lote. En particular, el porcentaje de conteos asignados a genes `protein_coding` tiene **mínimo 78,21 %** y **mediana 97,33 %** (archivo `02_estadisticos_QC.csv`): esa observación extrema requiere revisión individual y de posibles factores técnicos, **sin eliminar la muestra automáticamente**. Los conteos génicos sumados no equivalen al total de lecturas secuenciadas.

Evidencias: [`01_resumen_E05.csv`](docs/E_Analisis_Bulk/05_expresion/01_resumen_E05.csv), [`02_estadisticos_QC.csv`](docs/E_Analisis_Bulk/05_expresion/02_estadisticos_QC.csv), [`03_correspondencia_genes_GSE176078.csv`](docs/E_Analisis_Bulk/05_expresion/03_correspondencia_genes_GSE176078.csv), [`04_metodologia_E05.md`](docs/E_Analisis_Bulk/05_expresion/04_metodologia_E05.md) y [`05_sessionInfo.txt`](docs/E_Analisis_Bulk/05_expresion/05_sessionInfo.txt).

Figuras: [distribución de conteos](figures/E_Analisis_Bulk/05_expresion/01_conteos_por_muestra_hist.png) y [genes detectados](figures/E_Analisis_Bulk/05_expresion/02_genes_detectados_hist.png).

## 5. Entorno de trabajo y reproducibilidad

| Componente | Versión / implementación observada |
|---|---|
| Sistema operativo y editor | Windows 11, Positron |
| R | **4.5.3** (Windows UCRT x64) |
| Seurat | **5.6.0** |
| SeuratObject | **5.4.0** |
| BPCells | **0.3.1** |
| renv | **1.3.1** |
| Transferencia TCGA-BRCA | **GDC Data Transfer Tool 2.3** (PowerShell) |
| Versionado | Git y GitHub, rama `main` |
| Registro de paquetes | [`renv.lock`](renv.lock) |
| Evidencia de ejecución E/05 | [`05_sessionInfo.txt`](docs/E_Analisis_Bulk/05_expresion/05_sessionInfo.txt) |

**Estado de reproducibilidad:** los scripts E/01–E/05, los resúmenes agregados y los catálogos GDC públicos están versionados, y los manifiestos y MD5 permiten auditar los datos de entrada. Sin embargo, **aún no está probada** la reconstrucción íntegra desde un clon limpio mediante `renv::restore()`. Tampoco está demostrada la portabilidad de los `.rds` que dependen del almacén externo de BPCells. El flujo A requiere reconstruir dichos objetos desde los insumos originales si sus rutas locales no son válidas en otro equipo.

Los datos brutos (`data/raw/`) y las matrices, objetos Seurat/BPCells y métricas detalladas por muestra (`data/processed/`) **se excluyen de Git** con `.gitignore`. **Excepción deliberada:** las tablas de inventario y selección E/01–E/03 publicadas bajo `docs/E_Analisis_Bulk/` sí incluyen identificadores **públicos del GDC** (`Case ID`, `Sample ID` y/o `File ID`) para auditar el piloto. No se publican las matrices de expresión, TSV clínicos originales ni la tabla clínica paciente-nivel procesada. La ausencia de matrices en Git obliga a obtener de nuevo los archivos públicos originales para reproducir análisis desde cero; el informe puede regenerarse usando los CSV versionados. No se publican credenciales ni rutas personales.

## 6. Procedimiento de reproducción en Windows y Positron

> Todos los comandos deben ejecutarse **desde la raíz del repositorio**. `source()` corresponde a la **consola de R**; `git` y `gdc-client` corresponden a **PowerShell**. Los scripts R E/01–E/05 usan R base y utilidades incluidas en R.

### Paso 0. Obtener el repositorio y los insumos originales

```powershell
git clone https://github.com/5leyla5/proyecto-RNAseq-bulk.git
cd proyecto-RNAseq-bulk
```

Abrir la carpeta clonada en Positron. Obtener manualmente los archivos originales de SCP1039/GSE176078 y las exportaciones GDC de TCGA-BRCA descritas en la sección 2. Para repetir **exactamente el piloto del corte de noviembre de 2026**, usar los datos de entrada y el manifiesto versionados, no iniciar una nueva selección desde el portal.

Comprobar insumos en la **consola R**:

```r
stopifnot(
  file.exists("data/raw/GSE176078/Allcells_raw_count_out.zip"),
  file.exists("data/raw/GSE176078/Whole_miniatlas_meta.csv"),
  file.exists("data/raw/TCGA_BRCA/gdc_sample_sheet.tsv"),
  file.exists("data/raw/TCGA_BRCA/gdc_manifest.txt"),
  file.exists("data/raw/TCGA_BRCA/cohort_TCGA_BRCA_LSQ.2026-11-16.tsv"),
  file.exists("data/raw/TCGA_BRCA/clinical/clinical.tsv")
)
```

Verificar los MD5 de GSE176078 con [`docs/01_registro_datos.csv`](docs/01_registro_datos.csv). La ruta y el método de obtención de los TSV clínicos de TCGA-BRCA se detallan en el protocolo E/02. No utilizar TSV de otro proyecto TCGA.

### Paso 1. Preparar cohorte y clínica (R)

```r
source("scripts/E_Analisis_Bulk/01_preparar_cohorte_TCGA_BRCA.R")
source("scripts/E_Analisis_Bulk/02_validar_clinica_TCGA_BRCA.R")
```

E/01 genera catálogo completo, selección de una muestra por paciente y piloto aleatorio reproducible. E/02 cruza clínica y consolidación por paciente. Sus resultados están en `docs/E_Analisis_Bulk/01_cohorte/` y `02_clinica/`. Los archivos por paciente quedan locales bajo `data/processed/TCGA_BRCA/`.

### Paso 2. Cerrar la selección del piloto (R)

```r
source("scripts/E_Analisis_Bulk/03_cerrar_piloto_TCGA_BRCA.R")
```

E/03 preserva exactamente los 50 identificadores seleccionados en E/01, comprueba los estadios y genera [`03_manifest_piloto_50_FINAL.txt`](docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt).

### Paso 3. Descargar la expresión del piloto (PowerShell)

Instalar/descomprimir `gdc-client` para Windows desde [NCI GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool). El ejecutable puede no estar en `PATH` y necesita una carpeta de destino existente. El método **efectivamente usado** en Windows fue:

```powershell
$gdc = Get-ChildItem "$env:USERPROFILE\Downloads" -Recurse -File -Filter "gdc-client*.exe" | Select-Object -First 1
if ($null -eq $gdc) { throw "No se encuentra gdc-client.exe. Indica la ruta de instalación." }
& $gdc.FullName --version
New-Item -ItemType Directory -Force "data/raw/TCGA_BRCA/expresion_piloto" | Out-Null
& $gdc.FullName download `
  -m "docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt" `
  -d "data/raw/TCGA_BRCA/expresion_piloto"
```

También existe un **script de automatización preparado, no probado en el equipo original**, que verifica el manifiesto, busca el ejecutable y registra su versión y hashes:

```powershell
& ".\scripts\E_Analisis_Bulk\04_descargar_bulk_TCGA_BRCA.ps1"
# Si no detecta el ejecutable:
# & ".\scripts\E_Analisis_Bulk\04_descargar_bulk_TCGA_BRCA.ps1" -GdcExe "C:\ruta\gdc-client.exe"
```

Repetir el mismo manifiesto y destino permite reanudar una descarga interrumpida. **No descargar el manifiesto global de 1.111 archivos** si el objetivo es reproducir el piloto de 50.

### Paso 4. Comprobar integridad de descarga (R)

```r
source("scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R")
```

Criterio de aceptación observado: **50 presentes, 50 con MD5 y tamaño correctos, 0 fallos**. La ejecución genera `docs/E_Analisis_Bulk/04_descarga/01_integridad_archivos.csv` y `02_resumen_descarga.csv`. Las matrices crudas permanecen en `data/raw/TCGA_BRCA/expresion_piloto/<UUID>/<archivo>.tsv`.

### Paso 5. Leer matrices y obtener QC preliminar (R)

```r
source("scripts/E_Analisis_Bulk/05_leer_QC_expresion_TCGA_BRCA.R")
```

El script E/05 requiere los 50 archivos, sus MD5 correctos y el ZIP de GSE176078. Guarda tablas agregadas en `docs/E_Analisis_Bulk/05_expresion/`, figuras en `figures/E_Analisis_Bulk/05_expresion/` y las matrices locales bajo `data/processed/TCGA_BRCA/05_expresion/` (fuera de Git). Una ejecución que no encuentre la referencia informa el problema: **no** debe interpretarse un solapamiento faltante como cero genes comunes.

### Paso 6. Reproducir la referencia scRNA (objetivo A)

Revisar los scripts disponibles bajo `scripts/A_Analisis_SingleCell/` y los registros del directorio `docs/`. Los insumos originales son GSE176078/SCP1039. El QC inicial A se ejecutó en el equipo de origen; sus objetos `.rds` y el almacén BPCells son locales. **La restauración completa de A desde un clon limpio aún está pendiente de prueba y documentación detallada.**

## 7. Flujo de análisis A–F y estado actual

| Objetivo | Trabajo | Estado al 17-11-2026 |
|---|---|---|
| **A. scRNA-seq** | Acceso, carga, controles y QC inicial GSE176078 | **Ejecutado:** QC descriptivo inicial |
| **B. Referencia celular** | Separar donantes antes de selección de genes y construir perfiles celulares | **Pendiente** |
| **C. Pseudobulk y deconvolución** | Crear mezclas de donantes reservados; NNLS/MuSiC y comparador basal | **Pendiente** |
| **D. Evaluación** | MAE, RMSE y errores por tipo celular, evitando pseudorreplicación | **Pendiente** |
| **E/01–E/03. Cohorte TCGA** | Selección, clínica y manifiesto de 50 pacientes | **Ejecutado** |
| **E/04. Descarga y MD5** | Obtener y auditar los 50 archivos STAR Counts | **Ejecutado y verificado** |
| **E/05. QC bulk** | Lectura de 50 matrices, QC descriptivo y solapamiento preliminar de símbolos | **Ejecutado** |
| **E posterior. Aplicación** | Armonizar genes, aplicar referencia validada a TCGA-BRCA | **Pendiente** |
| **F. Discusión** | Sesgos, validez, alcance y limitaciones | **Pendiente** |

**Diseño de validación propuesto, todavía no ejecutado:** partición aproximada **70/30 por donante** (18 entrenamiento y 8 prueba, sujeta a distribución de ER+/HER2+/TNBC y cobertura celular), con semilla y asignación documentadas. La selección de genes y construcción de firmas se hará **solo sobre entrenamiento**. Las mezclas pseudobulk de prueba sumarán conteos UMI crudos de células de donantes reservados; se conservarán la composición celular conocida, el total de células y la identidad del donante. Se compararán NNLS y, si es factible, MuSiC contra un basal de composición promedio del entrenamiento, con MAE, RMSE y sesgos por tipo. La unidad independiente de resumen será **el donante**: varias mezclas de una misma persona no son réplicas biológicas independientes. No se utilizarán los donantes de evaluación para optimizar marcadores.

**Precauciones:** TCGA-BRCA no aporta proporciones celulares verdaderas para validar las predicciones. La asociación de composición celular estimada con estadio AJCC será observacional y **exploratoria**, incluso si posteriormente se amplía la cohorte; los 20 pacientes estadio IV de la cohorte completa y la categoría de preservación `Unknown` limitan la interpretación. No confundir asociación con causalidad.

## 8. Informe R Markdown y evidencias de la Entrega 1

Los **cinco criterios de diseño y factibilidad** se respaldan así:

| Criterio | Evidencia auditada | Situación |
|---|---|---|
| Selección y acceso a datos | GSE176078 y GDC TCGA-BRCA; MD5, catálogos y 50 matrices | **Demostrado** |
| Definición de cohorte | E/01, E/02 y E/03; criterios de selección y 50 casos finales | **Demostrado** |
| Unidades de análisis | Donante/célula en scRNA; paciente/muestra en bulk | **Definidas** |
| Tabla preliminar de muestras | Resúmenes y QC E/01–E/05; tabla de sección 3 | **Demostrado** |
| Ambiente computacional inicial | R/Positron, `renv.lock`, `sessionInfo`, herramienta GDC | **Configurado; restauración limpia pendiente** |

**Informe reproducible ya publicado:**

- [Fuente R Markdown — `Informe_Entrega_1.Rmd`](docs/Entrega_1/Informe_Entrega_1.Rmd): texto, verificaciones automáticas y cifras obtenidas de los CSV versionados.
- [Informe HTML generado — `Informe_Entrega_1.html`](docs/Entrega_1/Informe_Entrega_1.html): salida de lectura/impresión; incluye la imagen del flujo incrustada para que no dependa de rutas locales.
- [Diagrama del flujo](docs/Entrega_1/assets/flujo_metodologico.png) y [estilos de impresión](docs/Entrega_1/estilo_informe.css).
- [Script de generación](scripts/00_render_informe_entrega1.R) y [procedimiento de reproducción](docs/Entrega_1/LEEME_informe.md).

**Regeneración (consola R de Positron, desde la raíz):**

```r
source("scripts/00_render_informe_entrega1.R")
```

El informe se regenera a partir de **resultados tabulados ya versionados**, sin necesitar descargar las 50 matrices cada vez. No equivale a haber probado la reconstrucción íntegra del análisis original en un equipo nuevo. Antes de entregar, comprobar en la vista de impresión que el HTML quede en **tres páginas A4**, con encabezado/pie del navegador desactivados si corresponde, y verificar que imagen, tablas, resultados y bibliografía se vean íntegros. No se ha generado un Word; no es necesario ejecutar MuSiC, NNLS ni inferencias clínicas para esta etapa.

## 9. Limitaciones y tareas futuras

- **Armonización de genes:** resolver Ensembl con versión frente a símbolos de GSE176078, alias, duplicados y anotaciones. El valor **20.288** es una comprobación preliminar, no una matriz armonizada.
- **QC:** las distribuciones E/05 no bastan para declarar todas las muestras de calidad equivalente. Investigar la muestra con 78,21 % de conteos `protein_coding` (mediana 97,33 %), posibles lotes y predefinir criterios de filtrado antes de excluir casos.
- **Entrenamiento y validación:** verificar cobertura desigual de categorías por donante (13/26 Normal Epithelial; 14/26 Plasmablasts), mantener a los donantes de evaluación fuera de selección de genes y ajuste, y resumir errores por donante, no por célula o mezcla.
- **Entorno:** probar `renv::restore()` desde un clon limpio antes de declarar reproducibilidad integral; revisar dependencias de BPCells y rutas locales.
- **Análisis clínico:** el piloto no está diseñado para comparaciones por estadio; contemplar ampliar a la cohorte no FFPE tras el QC y comprobar cuántos pacientes tienen simultáneamente expresión elegible y estadio conocido.
- **Referencias metodológicas:** MuSiC y NNLS son métodos **propuestos**. No hay todavía proporciones estimadas, métricas de desempeño ni asociaciones clínicas.

## 10. Referencias y documentación

1. Wu SZ et al. (2021). *A single-cell and spatially resolved atlas of human breast cancers*. **Nature Genetics**. https://doi.org/10.1038/s41588-021-00911-1
2. Wang X et al. (2019). *Bulk tissue cell type deconvolution with multi-subject single-cell expression reference*. **Nature Communications** 10:380. https://doi.org/10.1038/s41467-018-08023-x
3. Avila Cobos F et al. (2020). *Benchmarking of cell type deconvolution pipelines for transcriptomics data*. **Nature Communications**. https://doi.org/10.1038/s41467-020-19015-1
4. Tran KA et al. (2023). *Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures*. **Nature Communications**. https://doi.org/10.1038/s41467-023-41385-5
5. Luecken MD & Theis FJ (2019). *Current best practices in single-cell RNA-seq analysis: a tutorial*. **Molecular Systems Biology**. https://doi.org/10.15252/msb.20188746
6. [MuSiC: documentación/tutorial oficial para R](https://xuranw.github.io/MuSiC/articles/MuSiC.html).
7. [NCI GDC — Guías oficiales](https://docs.gdc.cancer.gov/), especialmente [STAR Counts/mRNA pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) y [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/).

---

**Integridad científica:** las cifras documentadas provienen de CSV generados por scripts ejecutados y versionados. B–D y la aplicación final de E siguen pendientes. No se distribuyen matrices crudas ni datos clínicos originales por paciente en GitHub; **sí** se distribuyen tablas de metadatos con identificadores públicos GDC necesarios para auditar la selección, como se aclara en la sección 5.
