# =============================================================================
# BIOINFORMATICA II | Proyecto RNA-seq bulk | Objetivo E / etapa 01
# Catalogo reproducible de TCGA-BRCA y seleccion piloto (sin datos clinicos aun)
# Fuente: archivos exportados desde GDC Repository, 2026-11-16.
# Ejecutar en la CONSOLE de R, desde la raiz del proyecto:
#   source("scripts/E_Analisis_Bulk/01_preparar_cohorte_TCGA_BRCA.R")
# No descarga expresion, no requiere paquetes adicionales y no altera GSE176078.
# =============================================================================

if (!file.exists("README.md") || !dir.exists("scripts/E_Analisis_Bulk")) {
  stop("Abre en Positron la carpeta raiz 'Proyecto RNA-seq bulk'.")
}

origen <- file.path("data", "raw", "TCGA_BRCA")
sheet_path <- file.path(origen, "gdc_sample_sheet.tsv")
manifest_path <- file.path(origen, "gdc_manifest.txt")
tabla_path <- file.path(origen, "files-table.tsv")
salidas <- file.path("docs", "E_Analisis_Bulk", "01_cohorte")
dir.create(salidas, recursive = TRUE, showWarnings = FALSE)

faltantes <- c(sheet_path, manifest_path)
if (!all(file.exists(faltantes))) {
  stop("Archivos GDC faltantes: ", paste(faltantes[!file.exists(faltantes)],
                                       collapse = ", "),
       ". Ver README: descargar Sample Sheet y Manifest de GDC Repository.")
}

leer_tsv <- function(archivo) {
  utils::read.delim(archivo, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8", quote = "")
}
s <- leer_tsv(sheet_path)
m <- leer_tsv(manifest_path)

campos_s <- c("File ID", "File Name", "Project ID", "Data Type",
              "Case ID", "Sample ID", "Tissue Type", "Tumor Descriptor",
              "Preservation Method")
campos_m <- c("id", "filename", "md5", "size", "state")
if (!all(campos_s %in% names(s))) {
  stop("Columnas ausentes en Sample Sheet: ",
       paste(setdiff(campos_s, names(s)), collapse = ", "))
}
if (!all(campos_m %in% names(m))) {
  stop("Columnas ausentes en Manifest: ",
       paste(setdiff(campos_m, names(m)), collapse = ", "))
}
if (anyDuplicated(s[["File ID"]]) || anyDuplicated(m$id)) {
  stop("Hay identificadores de archivos duplicados: revisar exportacion GDC.")
}
if (!setequal(s[["File ID"]], m$id)) {
  stop("Sample Sheet y Manifest no tienen los mismos File ID. No continuar.")
}
ind <- match(s[["File ID"]], m$id)
if (!identical(s[["File Name"]], m$filename[ind])) {
  stop("Los nombres de archivo difieren entre Sample Sheet y Manifest.")
}
if (!all(s[["Project ID"]] == "TCGA-BRCA")) {
  stop("El Sample Sheet contiene proyectos distintos de TCGA-BRCA.")
}
if (!all(s[["Data Type"]] == "Gene Expression Quantification")) {
  stop("El Sample Sheet contiene datos distintos de expresion genica.")
}
if (!all(tolower(s[["Tissue Type"]]) == "tumor") ||
    !all(tolower(s[["Tumor Descriptor"]]) == "primary")) {
  stop("La exportacion incluye muestras distintas de tumor primario.")
}
if (!all(substr(s[["Sample ID"]], 14L, 15L) == "01")) {
  stop("Algun barcode no indica tumor primario (codigo TCGA 01). Revisar.")
}
if (!all(substr(s[["Sample ID"]], 1L, 12L) == s[["Case ID"]])) {
  stop("Hay discordancia entre identificador de paciente y muestra.")
}

# Unir los checksums y los tamanos SIN depender del orden de exportacion.
s$md5 <- m$md5[ind]
s$size_bytes <- as.numeric(m$size[ind])
s$state <- m$state[ind]
s <- s[order(s[["Case ID"]], s[["Sample ID"]], s[["File ID"]]),
       , drop = FALSE]
rownames(s) <- NULL

# Seleccion de UN archivo por paciente.
# Regla preliminar y determinista, NO sustituye una revision tecnica de QC:
# 1) favorecer preservacion conocida OCT; 2) Unknown; 3) FFPE;
# desempatar por Sample ID y File ID alfabeticos.
s$prioridad_preservacion <- match(s[["Preservation Method"]],
                                  c("OCT", "Unknown", "FFPE"))
s$prioridad_preservacion[is.na(s$prioridad_preservacion)] <- 99L
s_orden <- s[order(s[["Case ID"]], s$prioridad_preservacion,
                    s[["Sample ID"]], s[["File ID"]]), , drop = FALSE]
una_por_paciente <- s_orden[!duplicated(s_orden[["Case ID"]]), , drop = FALSE]
stopifnot(!anyDuplicated(una_por_paciente[["Case ID"]]))

# Para el PILOTO excluir FFPE: reduce una fuente potencial de efecto de lote.
# Guardar la cohorte completa por separado para no esconder exclusiones.
candidatos <- una_por_paciente[
  una_por_paciente[["Preservation Method"]] != "FFPE", , drop = FALSE
]
if (nrow(candidatos) < 50L) stop("Hay menos de 50 pacientes no FFPE.")

# Seleccion aleatoria reproducible, TODAVIA PROVISIONAL: faltan datos clinicos.
set.seed(20261116L)
seleccion <- sort(sample.int(nrow(candidatos), 50L, replace = FALSE))
piloto <- candidatos[seleccion, , drop = FALSE]
piloto <- piloto[order(piloto[["Case ID"]]), , drop = FALSE]
stopifnot(nrow(piloto) == 50L, !anyDuplicated(piloto[["Case ID"]]))

# Un manifiesto de GDC utiliza exactamente estos campos y este orden.
manifest_piloto <- data.frame(
  id = piloto[["File ID"]],
  filename = piloto[["File Name"]],
  md5 = piloto$md5,
  size = piloto$size_bytes,
  state = piloto$state,
  stringsAsFactors = FALSE
)

cuantos_por_caso <- table(s[["Case ID"]])
pacientes_multiples <- sum(cuantos_por_caso > 1L)
solo_ffpe <- sum(una_por_paciente[["Preservation Method"]] == "FFPE")
resumen <- data.frame(
  indicador = c(
    "Archivos STAR Counts exportados (tumor primario)",
    "Pacientes distintos (Case ID)",
    "Identificadores Sample ID distintos",
    "Pacientes con multiples archivos",
    "Archivos adicionales por duplicacion de pacientes",
    "Archivos con preservacion OCT",
    "Archivos con preservacion Unknown",
    "Archivos con preservacion FFPE",
    "Pacientes con solo archivos FFPE (segun regla de seleccion)",
    "Pacientes no FFPE elegibles para piloto",
    "Pacientes seleccionados en piloto provisional",
    "Tamano del catalogo completo (GB decimales)",
    "Filas visibles en files-table.tsv (no catalogo completo)"
  ),
  valor = c(
    nrow(s), length(unique(s[["Case ID"]])),
    length(unique(s[["Sample ID"]])), pacientes_multiples,
    nrow(s) - length(unique(s[["Case ID"]])),
    sum(s[["Preservation Method"]] == "OCT"),
    sum(s[["Preservation Method"]] == "Unknown"),
    sum(s[["Preservation Method"]] == "FFPE"),
    solo_ffpe, nrow(candidatos), nrow(piloto),
    round(sum(s$size_bytes) / 1e9, 3),
    if (file.exists(tabla_path)) nrow(leer_tsv(tabla_path)) else NA_real_
  ),
  fuente = "GDC Repository exportado el 2026-11-16",
  stringsAsFactors = FALSE
)

utils::write.csv(resumen, file.path(salidas, "01_resumen_TCGA_BRCA.csv"),
                 row.names = FALSE)
utils::write.csv(s, file.path(salidas, "02_catalogo_TCGA_BRCA.csv"),
                 row.names = FALSE)
utils::write.csv(una_por_paciente,
                 file.path(salidas, "03_una_muestra_por_paciente.csv"),
                 row.names = FALSE)
utils::write.csv(piloto, file.path(salidas, "04_piloto_50_provisional.csv"),
                 row.names = FALSE)
utils::write.table(manifest_piloto,
                   file.path(salidas, "05_manifest_piloto_50.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE,
                   col.names = TRUE)

notas <- c(
  "# TCGA-BRCA: definicion preliminar de cohorte",
  "",
  "Fuente: GDC Data Portal, Repository (consulta/descarga 2026-11-16).",
  "Proyecto TCGA-BRCA; Gene Expression Quantification; RNA-Seq; STAR - Counts;",
  "Transcriptome Profiling; acceso Open; Tissue Type Tumor; Tumor Descriptor Primary.",
  "Se comprobo adicionalmente codigo TCGA 01 en cada Sample ID.",
  "",
  "## Unidad de analisis y elegibilidad",
  "Una unidad clinica = un Case ID; una muestra de expresion por paciente.",
  "En pacientes con varios archivos se priorizo OCT, luego Unknown, luego FFPE;",
  "los empates se resolvieron por Sample ID y File ID. Es una regla preliminar",
  "que debera revisarse tras examinar la calidad y el efecto de lote.",
  "",
  "El piloto provisional de 50 pacientes excluye FFPE y utiliza set.seed(20261116).",
  "La seleccion debera revisarse tras obtener datos clinicos (estadio y faltantes).",
  "Las 1111 matrices de expresion NO han sido descargadas por este script.",
  "Los metadatos originales permanecen en data/raw/TCGA_BRCA/ (excluido de Git).",
  "",
  "## Comprobaciones reproducibles",
  "Se verifico igualdad por File ID y nombre entre Sample Sheet y Manifest;",
  "tambien se verificaron proyecto, tipo de dato, tipo y descriptor tumoral,",
  "codigo de muestra 01 y coincidencia Case ID / Sample ID.",
  "",
  "## Pendiente",
  "Descargar clinical.tsv/clinical.json del GDC, evaluar estadio AJCC y faltantes,",
  "confirmar seleccion piloto, descargar y verificar MD5 de expresion,",
  "comprobar correspondencia de genes con GSE176078 y hacer deconvolucion."
)
writeLines(notas, file.path(salidas, "06_notas_metodologicas.md"), useBytes = TRUE)

cat("\n== TCGA-BRCA / E01: CATALOGO VERIFICADO ==\n")
print(resumen[, c("indicador", "valor")], row.names = FALSE)
cat("\nSalidas reproducibles en: ", salidas, "\n", sep = "")
