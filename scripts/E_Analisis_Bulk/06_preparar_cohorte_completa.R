# Bioinformatica II - TCGA-BRCA
# E/06: Ampliacion a cohorte no FFPE completa
# No descarga ni modifica datos originales.

entrada <- "docs/E_Analisis_Bulk/01_cohorte/03_una_muestra_por_paciente.csv"
piloto <- "docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt"
salida <- "docs/E_Analisis_Bulk/06_cohorte_completa"

stopifnot(file.exists(entrada), file.exists(piloto))

datos <- read.csv(
  entrada, check.names = FALSE,
  stringsAsFactors = FALSE
)

requeridas <- c(
  "File ID", "File Name", "Case ID",
  "Sample ID", "Project ID",
  "Tissue Type", "Tumor Descriptor",
  "Preservation Method", "md5",
  "size_bytes", "state"
)

if (!all(requeridas %in% names(datos))) {
  stop("Faltan columnas en el catalogo")
}

stopifnot(
  nrow(datos) == 1095L,
  !anyDuplicated(datos[["Case ID"]]),
  all(datos[["Project ID"]] == "TCGA-BRCA"),
  all(datos[["Tissue Type"]] == "Tumor"),
  all(datos[["Tumor Descriptor"]] == "Primary"),
  all(grepl("-01", datos[["Sample ID"]], fixed = TRUE)),
  all(datos[["Preservation Method"]] %in%
        c("OCT", "Unknown", "FFPE"))
)

cohorte <- datos[
  datos[["Preservation Method"]] %in% c("OCT", "Unknown"),
  , drop = FALSE
]

stopifnot(
  nrow(cohorte) == 1093L,
  !anyDuplicated(cohorte[["Case ID"]]),
  !anyDuplicated(cohorte[["File ID"]]),
  all(cohorte[["state"]] == "released"),
  all(grepl("^[0-9a-fA-F]{32}$", cohorte[["md5"]])),
  all(is.finite(cohorte[["size_bytes"]])),
  all(cohorte[["size_bytes"]] > 0)
)

manifest <- data.frame(
  id = cohorte[["File ID"]],
  filename = cohorte[["File Name"]],
  md5 = cohorte[["md5"]],
  size = cohorte[["size_bytes"]],
  state = cohorte[["state"]]
)

p <- read.delim(piloto, stringsAsFactors = FALSE)
stopifnot(
  nrow(p) == 50L,
  all(p$id %in% manifest$id)
)

dir.create(salida, recursive = TRUE, showWarnings = FALSE)

write.table(
  manifest,
  file.path(salida, "manifest_1093_TCGA_BRCA.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.csv(
  cohorte,
  file.path(salida, "cohorte_1093_seleccionada.csv"),
  row.names = FALSE
)

resumen <- data.frame(
  indicador = c(
    "Pacientes catalogo",
    "Pacientes no FFPE",
    "Archivos OCT",
    "Archivos Unknown",
    "Piloto original incluido",
    "Tamano total GB"
  ),
  valor = c(
    nrow(datos),
    nrow(cohorte),
    sum(cohorte[["Preservation Method"]] == "OCT"),
    sum(cohorte[["Preservation Method"]] == "Unknown"),
    nrow(p),
    round(sum(manifest$size) / 1e9, 3)
  )
)

write.csv(
  resumen,
  file.path(salida, "resumen_cohorte_completa.csv"),
  row.names = FALSE
)

cat("\n=== COHORTE COMPLETA TCGA-BRCA ===\n")
print(resumen, row.names = FALSE)
cat("\nManifest validado: ", nrow(manifest), " archivos\n")
