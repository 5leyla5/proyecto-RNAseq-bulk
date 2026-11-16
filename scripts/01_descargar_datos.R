# ============================================================
# BIOINFORMATICA II - Proyecto RNA-seq bulk
# 01 - Registrar y verificar datos originales (GSE176078)
# Fuente: Broad Single Cell Portal, estudio SCP1039
# Ejecutar desde la raiz: source("scripts/01_descargar_datos.R")
# NO descarga otra copia, NO descomprime ni filtra los datos.
# ============================================================

# 1. Rutas relativas: funcionan en Windows, Linux y macOS.
if (!file.exists("README.md") || !dir.exists("scripts")) {
  stop("Abre la carpeta raiz 'Proyecto RNA-seq bulk' en Positron.")
}

carpeta <- file.path("data", "raw", "GSE176078")
ruta_zip <- file.path(carpeta, "Allcells_raw_count_out.zip")
ruta_meta <- file.path(carpeta, "Whole_miniatlas_meta.csv")
archivos <- c(ruta_zip, ruta_meta)

if (!all(file.exists(archivos))) {
  stop("Faltan originales en data/raw/GSE176078/: ",
       paste(basename(archivos[!file.exists(archivos)]), collapse = ", "),
       "\nDescargalos desde https://singlecell.broadinstitute.org/single_cell/study/SCP1039")
}

# 2. Verificar que el ZIP contiene exactamente el formato 10x esperado.
inventario <- utils::unzip(ruta_zip, list = TRUE)
componentes <- paste0("BrCa_Atlas_Count_out/",
                     c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
stopifnot(all(componentes %in% inventario$Name))

# 3. Validar identificadores de celulas presentes en metadatos.
metadatos <- utils::read.csv(ruta_meta, stringsAsFactors = FALSE,
                             check.names = FALSE)
stopifnot(all(c("NAME", "Patient", "celltype_major") %in% names(metadatos)))
# Broad incluye a veces una fila TYPE con tipos de columnas, no es una celula.
if (nrow(metadatos) > 0L && metadatos$NAME[1] == "TYPE") {
  metadatos <- metadatos[-1L, , drop = FALSE]
}
stopifnot(nrow(metadatos) == 100064L,
          !anyNA(metadatos$NAME),
          anyDuplicated(metadatos$NAME) == 0L)

# 4. Guardar procedencia e integridad para reproducibilidad.
dir.create("docs", recursive = TRUE, showWarnings = FALSE)
fuente <- "https://singlecell.broadinstitute.org/single_cell/study/SCP1039"
registro <- data.frame(
  estudio = "GSE176078 / SCP1039",
  archivo = basename(archivos),
  origen = fuente,
  bytes = unname(file.info(archivos)$size),
  md5 = unname(tools::md5sum(archivos)),
  stringsAsFactors = FALSE
)

# Si hubo una verificacion anterior, avisar si cambiaron los originales.
ruta_registro <- file.path("docs", "01_registro_datos.csv")
if (file.exists(ruta_registro)) {
  anterior <- utils::read.csv(ruta_registro, stringsAsFactors = FALSE)
  if (all(c("archivo", "md5") %in% names(anterior))) {
    anterior_md5 <- setNames(anterior$md5, anterior$archivo)[registro$archivo]
    if (any(!is.na(anterior_md5) & anterior_md5 != registro$md5)) {
      stop("Cambio el MD5 de un archivo original: revisar antes de continuar.")
    }
  }
}
utils::write.csv(registro, ruta_registro, row.names = FALSE)

# 5. Confirmar que Git no suba archivos grandes por accidente.
ruta_ignore <- ".gitignore"
lineas <- if (file.exists(ruta_ignore)) readLines(ruta_ignore, warn = FALSE) else character()
if (!"data/raw/" %in% lineas) {
  writeLines(c(lineas, "data/raw/"), ruta_ignore)
}

cat("\n== DATOS ORIGINALES GSE176078 VERIFICADOS ==\n")
print(registro, row.names = FALSE)
cat("Celulas en metadatos:", nrow(metadatos), "\n")
cat("Componentes 10x en el ZIP:", sum(componentes %in% inventario$Name), "de 3\n")
cat("Registro:", ruta_registro, "\n")
cat("PASO 01 COMPLETADO. No se modificaron los originales.\n")
