# =====================================================================
# BIOINFORMATICA II - Proyecto RNA-seq bulk / Objetivo A
# 04 - Documentar decision de QC y filtrado ADICIONAL
#
# Entradas:
#   data/processed/GSE176078/seurat_qc_inicial.rds
#   docs/A_Analisis_SingleCell/03_qc_por_grupos/
#       qc_genes_impacto_umbrales_pre.csv
#
# Salidas:
#   docs/A_Analisis_SingleCell/04_filtrado_qc/
#       decision_filtrado.md
#       resumen_conservacion.csv
#       umbrales_genes_evaluados.csv
#       celulas_por_tipo.csv
#       celulas_por_donante.csv
#       celulas_por_donante_y_tipo.csv
#       parametros_decision.csv
#   logs/A_Analisis_SingleCell/04_filtrado_qc/session_info.txt
#
# NO elimina celulas ni genes, NO normaliza, NO modifica el RDS.
# Ejecutar desde la raiz del proyecto:
# source("scripts/A_Analisis_SingleCell/04_documentar_decision_qc.R")
# =====================================================================

if (!file.exists("scripts/00_inicializar_proyecto.R")) {
  stop("Abre la raiz de 'Proyecto RNA-seq bulk' en Positron.")
}
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("Falta el paquete SeuratObject.")
}

base_datos <- Sys.getenv("RNASEQ_DATA_DIR", unset = "")
if (!nzchar(base_datos)) base_datos <- file.path("data", "processed")

entrada_obj <- file.path(base_datos, "GSE176078", "seurat_qc_inicial.rds")
entrada_genes <- file.path(
  "docs", "A_Analisis_SingleCell", "03_qc_por_grupos",
  "qc_genes_impacto_umbrales_pre.csv"
)
salida <- file.path("docs", "A_Analisis_SingleCell", "04_filtrado_qc")
salida_log <- file.path("logs", "A_Analisis_SingleCell", "04_filtrado_qc")

for (ruta in c(entrada_obj, entrada_genes)) {
  if (!file.exists(ruta)) stop("Falta archivo de entrada: ", ruta)
}
dir.create(salida, recursive = TRUE, showWarnings = FALSE)
dir.create(salida_log, recursive = TRUE, showWarnings = FALSE)

# 1) Comprobar datos de entrada sin cargar toda la matriz a RAM.
obj <- readRDS(entrada_obj)
meta <- obj[[]]
campos <- c("Patient", "celltype_major", "nFeature_RNA", "nCount_RNA", "percent.mt")
if (!all(campos %in% names(meta))) {
  stop("Faltan metadatos: ", paste(setdiff(campos, names(meta)), collapse = ", "))
}
stopifnot(
  ncol(obj) == nrow(meta), ncol(obj) == 100064L,
  nrow(obj) == 29733L,
  !anyNA(meta[, campos, drop = FALSE]),
  identical(rownames(meta), colnames(obj)),
  length(unique(meta$Patient)) == 26L
)

# 2) Verificar umbrales para CELULAS documentados por Wu et al. (2021).
#    Ya se aplicaron en la elaboracion original del atlas; NO se reaplican.
verificacion_autores <- c(
  genes_mayor_200 = all(meta$nFeature_RNA > 200),
  umi_mayor_250 = all(meta$nCount_RNA > 250),
  mito_menor_20 = all(meta$percent.mt < 20)
)
if (!all(verificacion_autores)) {
  stop("Hay celulas fuera de los umbrales originales; revisar el estudio.")
}

# 3) Leer resultados de A/03, sin recalcular la matriz completa.
umbrales <- utils::read.csv(entrada_genes, stringsAsFactors = FALSE)
requeridas <- c(
  "umbral_min_celulas", "genes_conservados",
  "genes_excluidos_hipoteticos", "porcentaje_excluido_hipotetico"
)
if (!all(requeridas %in% names(umbrales))) {
  stop("Columnas inesperadas en A/03: ",
       paste(setdiff(requeridas, names(umbrales)), collapse = ", "))
}
umbrales <- umbrales[order(umbrales$umbral_min_celulas), requeridas]
if (!identical(as.integer(umbrales$umbral_min_celulas), c(1L, 3L, 10L))) {
  stop("Los escenarios de genes deben ser >=1, >=3 y >=10 celulas.")
}
if (!all(umbrales$genes_conservados +
         umbrales$genes_excluidos_hipoteticos == nrow(obj))) {
  stop("Los totales de escenarios de genes no concuerdan con el objeto.")
}
if (!all(abs(
    100 * umbrales$genes_excluidos_hipoteticos / nrow(obj) -
      umbrales$porcentaje_excluido_hipotetico
  ) < 0.02)) {
  stop("Porcentajes de exclusion incompatibles con el recuento de genes.")
}
if (umbrales$genes_excluidos_hipoteticos[1] != 0L) {
  stop("El escenario >=1 celula deberia conservar los genes existentes.")
}

# 4) Tablas de composicion ANTES de cualquier filtro adicional.
por_tipo <- as.data.frame(table(
  celltype_major = as.character(meta$celltype_major)),
  stringsAsFactors = FALSE
)
names(por_tipo)[2] <- "n_celulas_conservadas"
por_donante <- as.data.frame(table(
  Patient = as.character(meta$Patient)), stringsAsFactors = FALSE
)
names(por_donante)[2] <- "n_celulas_conservadas"
por_ambos <- as.data.frame(table(
  Patient = as.character(meta$Patient),
  celltype_major = as.character(meta$celltype_major)
), stringsAsFactors = FALSE)
names(por_ambos)[3] <- "n_celulas_conservadas"
stopifnot(
  sum(por_tipo$n_celulas_conservadas) == ncol(obj),
  sum(por_donante$n_celulas_conservadas) == ncol(obj),
  sum(por_ambos$n_celulas_conservadas) == ncol(obj)
)

resumen <- data.frame(
  indicador = c("celulas", "genes", "donantes", "tipos_celulares"),
  antes = c(ncol(obj), nrow(obj), nrow(por_donante), nrow(por_tipo)),
  eliminados_adicionales = 0L,
  despues = c(ncol(obj), nrow(obj), nrow(por_donante), nrow(por_tipo))
)
parametros <- data.frame(
  parametro = c(
    "decision", "filtro_adicional_celulas", "filtro_adicional_genes",
    "criterios_publicados_celulas", "escenarios_genes_solo_simulados",
    "unidad_metrica_mitocondrial", "objeto_entrada",
    "seleccion_marcadores_referencia"
  ),
  valor = c(
    "conservar_todos_los_registros_del_atlas_recibido", "ninguno", "ninguno",
    "nFeature_RNA>200; nCount_RNA>250; percent.mt<20",
    ">=1; >=3; >=10 celulas", "porcentaje_UMI_mitocondrial_por_celula",
    entrada_obj, "solo_tras_separar_donantes_objetivo_B"
  )
)

utils::write.csv(resumen, file.path(salida, "resumen_conservacion.csv"), row.names = FALSE)
utils::write.csv(umbrales, file.path(salida, "umbrales_genes_evaluados.csv"), row.names = FALSE)
utils::write.csv(por_tipo, file.path(salida, "celulas_por_tipo.csv"), row.names = FALSE)
utils::write.csv(por_donante, file.path(salida, "celulas_por_donante.csv"), row.names = FALSE)
utils::write.csv(por_ambos, file.path(salida, "celulas_por_donante_y_tipo.csv"), row.names = FALSE)
utils::write.csv(parametros, file.path(salida, "parametros_decision.csv"), row.names = FALSE)

# 5) Generar el documento metodologico (Markdown) a partir de los resultados.
formatear <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE,
                                trim = TRUE)
linea_gen <- function(i) {
  paste0("| >=", umbrales$umbral_min_celulas[i], " celulas | ",
         formatear(umbrales$genes_conservados[i]), " | ",
         formatear(umbrales$genes_excluidos_hipoteticos[i]), " | ",
         formatC(umbrales$porcentaje_excluido_hipotetico[i], format = "f",
                 digits = 2, decimal.mark = ","), " % |")
}
linea_tipo <- function(i) {
  paste0("| ", por_tipo$celltype_major[i], " | ",
         formatear(por_tipo$n_celulas_conservadas[i]), " |")
}

informe <- c(
  "# Objetivo A / 04 - Decision documentada de control de calidad",
  "",
  "## Decision",
  "Se conservaron las celulas y los genes del archivo scRNA-seq GSE176078 recibido",
  "del Broad Single Cell Portal, sin aplicar filtros adicionales en esta etapa.",
  "Los conteos crudos y los objetos previos NO fueron modificados.",
  "Esta es una decision fundamentada de QC, no la afirmacion de que todas",
  "las celulas son perfectas ni de que se descarta revisar posibles dobletes.",
  "",
  "## Punto de partida",
  paste0("- Celulas conservadas: **", formatear(ncol(obj)), "** (eliminadas adicionales: **0**)."),
  paste0("- Genes conservados: **", formatear(nrow(obj)), "** (eliminados adicionales: **0**)."),
  paste0("- Donantes: **", nrow(por_donante), "**."),
  paste0("- Tipos celulares principales: **", nrow(por_tipo), "**."),
  "- Capa utilizada: conteos UMI crudos del ensayo RNA, respaldados por BPCells.",
  "",
  "## Fundamento del control de calidad celular",
  "Wu et al. (2021) procesaron los datos originales con Cell Ranger y EmptyDrops,",
  "e indicaron >200 genes detectados, >250 UMI y <20 % mitocondrial.",
  "El objeto recibido ya satisface estos tres criterios en todas sus celulas:",
  paste0("- nFeature_RNA > 200: ", verificacion_autores[1], "."),
  paste0("- nCount_RNA > 250: ", verificacion_autores[2], "."),
  paste0("- percent.mt < 20: ", verificacion_autores[3], "."),
  "En A/02 se recalculo el porcentaje mitocondrial desde los conteos crudos,",
  "contrastandolo con los metadatos; en A/03 se compararon los indicadores",
  "por tipo celular y por donante.",
  "Las distribuciones son heterogeneas; un unico umbral adicional muy estricto",
  "podria excluir de manera desproporcionada poblaciones epiteliales tumorales",
  "u otras poblaciones de interes para la referencia de deconvolucion.",
  "",
  "## Evaluacion exploratoria del QC de genes",
  "Se definio un gen como detectado en una celula si el conteo UMI es >= 1.",
  "Los escenarios siguientes fueron simulados en A/03, pero NO aplicados:",
  "",
  "| Deteccion minima | Genes conservados | Genes excluidos hipoteticos | Porcentaje excluido |",
  "|---|---:|---:|---:|",
  vapply(seq_len(nrow(umbrales)), linea_gen, character(1)),
  "",
  "Se evaluaron ademas la distribucion de frecuencia de deteccion, la relacion",
  "entre UMI medio crudo y deteccion, y los 20 genes con mayor UMI acumulado.",
  "La baja deteccion global no basta para descartar un gen: podria ser relevante",
  "en una poblacion celular poco abundante.",
  "La seleccion de genes altamente variables (HVG) para PCA sera una etapa",
  "posterior y diferente del filtrado de genes del objeto original.",
  "",
  "## Celulas conservadas por tipo celular",
  "| Tipo celular | Celulas conservadas |",
  "|---|---:|",
  vapply(seq_len(nrow(por_tipo)), linea_tipo, character(1)),
  "",
  "Los recuentos por donante y por combinacion donante-tipo celular estan",
  "documentados en `celulas_por_donante.csv` y `celulas_por_donante_y_tipo.csv`.",
  "Incluyen las combinaciones sin celulas y no representan proporciones reales",
  "de tejido tumoral, sino distribucion de celulas capturadas y anotadas.",
  "",
  "## Limitaciones y decisiones futuras",
  "- El atlas ya tiene un preprocesamiento previo; no estamos evaluando gotas vacias originales.",
  "- El QC descriptivo no descarta por si solo dobletes, contaminacion ambiental o efectos de lote.",
  "- No se seleccionaran marcadores de referencia con datos de donantes de evaluacion.",
  "- En el objetivo B se separaran donantes antes de seleccionar los marcadores de la referencia.",
  "- Para pseudobulk se conservaran y emplearan conteos crudos UMI, no la capa normalizada.",
  "- Si aparecieran problemas posteriores, cualquier nuevo filtro se documentara",
  "  con nuevos parametros, recuentos por donante y tipo celular y una nueva version del analisis.",
  "",
  "## Reproducibilidad",
  "- Ejecutar desde la raiz del proyecto:",
  "  `source(\"scripts/A_Analisis_SingleCell/04_documentar_decision_qc.R\")`.",
  paste0("- Objeto de entrada: `", entrada_obj, "`."),
  paste0("- Resumen de genes de entrada: `", entrada_genes, "`."),
  "- Este script solo lee resultados de A/02 y A/03 y escribe tablas y este documento.",
  "- Fuente de los archivos originales y sumas MD5: `docs/01_registro_datos.csv`.",
  "- Versiones de software: `logs/A_Analisis_SingleCell/04_filtrado_qc/session_info.txt`.",
  "",
  "## Referencias",
  "1. Wu SZ et al. (2021). A single-cell and spatially resolved atlas of human breast cancers.",
  "   Nature Genetics. https://doi.org/10.1038/s41588-021-00911-1",
  "2. Single-cell best practices: Quality Control.",
  "   https://www.sc-best-practices.org/preprocessing-visualization/quality-control/",
  "",
  "**Decision registrada: no aplicar filtros adicionales a celulas ni genes.**"
)
writeLines(informe, file.path(salida, "decision_filtrado.md"), useBytes = TRUE)
writeLines(capture.output(utils::sessionInfo()), file.path(salida_log, "session_info.txt"))

cat("\n== OBJETIVO A/04: DECISION DE QC DOCUMENTADA ==\n")
print(resumen, row.names = FALSE)
cat("\nEscenarios de deteccion por gen:\n")
print(umbrales, row.names = FALSE)
cat("\nArchivos guardados en: ", salida, "\n", sep = "")
cat("NO se modificaron datos, objetos Seurat ni matrices BPCells.\n")
cat("PASO A/04 COMPLETADO.\n")
