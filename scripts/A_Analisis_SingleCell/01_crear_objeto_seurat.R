# =============================================================================
# BIOINFORMATICA II | Proyecto RNA-seq bulk | GSE176078
# OBJETIVO A / ETAPA 01 - Matriz gen-celula y objeto Seurat inicial (SIN filtros)
# =============================================================================
# Ejecutar desde la carpeta raiz del proyecto, en una sesion nueva de R:
#   source("scripts/A_Analisis_SingleCell/01_crear_objeto_seurat.R")
#
# ENTRADAS (verificadas en el paso 01):
#   data/raw/GSE176078/Allcells_raw_count_out.zip
#   data/raw/GSE176078/Whole_miniatlas_meta.csv
#
# SALIDAS:
#   data/processed/GSE176078/bpcells_counts/   (matriz almacenada en disco)
#   data/processed/GSE176078/seurat_inicial.rds
#   docs/A_Analisis_SingleCell/01_objeto_seurat/resumen.csv
#   docs/A_Analisis_SingleCell/01_objeto_seurat/versiones.csv
#   logs/A_Analisis_SingleCell/01_objeto_seurat/session_info.txt
#
# OPCIONAL: configurar RNASEQ_DATA_DIR con otra carpeta para salidas grandes.
# NO se seleccionan genes, NO se filtran celulas, NO se normaliza.
# =============================================================================

# 01A. Comprobar entradas y paquetes -----------------------------------------
ruta_zip  <- file.path("data", "raw", "GSE176078", "Allcells_raw_count_out.zip")
ruta_meta <- file.path("data", "raw", "GSE176078", "Whole_miniatlas_meta.csv")
ruta_reg  <- file.path("docs", "01_registro_datos.csv")

if (!all(file.exists(c(ruta_zip, ruta_meta, ruta_reg)))) {
  stop("Faltan datos originales o su registro. Ejecuta primero scripts/01_descargar_datos.R ",
       "desde la carpeta raiz del proyecto.")
}

paquetes <- c("BPCells", "Seurat", "SeuratObject", "Matrix")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Paquetes faltantes: ", paste(faltan, collapse = ", "))

# El valor por defecto es relativo al proyecto: funciona en Windows, Linux y macOS.
# Si no deseas archivos grandes en una carpeta sincronizada, configura RNASEQ_DATA_DIR.
base_datos <- Sys.getenv("RNASEQ_DATA_DIR", unset = "")
if (!nzchar(base_datos)) base_datos <- file.path("data", "processed")
base_datos <- file.path(base_datos, "GSE176078")

ruta_bp  <- file.path(base_datos, "bpcells_counts")
ruta_rds <- file.path(base_datos, "seurat_inicial.rds")
ruta_tmp <- file.path(base_datos, "_temporal_importacion")
ruta_10x <- file.path(ruta_tmp, "BrCa_Atlas_Count_out")

dir.create(base_datos, recursive = TRUE, showWarnings = FALSE)
dir.create("docs/A_Analisis_SingleCell/01_objeto_seurat", recursive = TRUE, showWarnings = FALSE)
dir.create("logs/A_Analisis_SingleCell/01_objeto_seurat", recursive = TRUE, showWarnings = FALSE)

if (grepl("onedrive", normalizePath(base_datos, winslash = "/", mustWork = FALSE),
          ignore.case = TRUE)) {
  warning("Los datos procesados estan dentro de OneDrive. Para mejorar el rendimiento, ",
          "puedes configurar RNASEQ_DATA_DIR en un disco local no sincronizado.")
}

# 01B. Crear o reutilizar la matriz BPCells ----------------------------------
if (!dir.exists(ruta_bp)) {
  componentes <- paste0(
    "BrCa_Atlas_Count_out/",
    c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  )
  inventario <- utils::unzip(ruta_zip, list = TRUE)
  if (!all(componentes %in% inventario$Name)) {
    stop("El ZIP no contiene la matriz 10x esperada.")
  }

  cat("\n[A/01] Extrayendo temporalmente los tres archivos 10x...\n")
  dir.create(ruta_tmp, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(ruta_zip, files = componentes, exdir = ruta_tmp, overwrite = TRUE)
  if (!all(file.exists(file.path(ruta_tmp, componentes)))) {
    stop("La extraccion no termino correctamente.")
  }

  cat("[A/01] Convirtiendo a BPCells con memoria limitada...\n")
  # Los conteos quedan en disco; evita crear una matriz densa en RAM.
  BPCells::import_matrix_market_10x(
    mtx_dir = ruta_10x,
    outdir = ruta_bp,
    tmpdir = ruta_tmp,
    sort_bytes = 268435456
  )
  # Borrar solo la copia TEMPORAL que acabamos de extraer; conservar originales.
  unlink(ruta_tmp, recursive = TRUE)
} else {
  cat("\n[A/01] Matriz BPCells existente: no se repite la importacion.\n")
}

matriz <- tryCatch(
  BPCells::open_matrix_dir(ruta_bp),
  error = function(e) stop("No se pudo abrir BPCells. Revisa si la carpeta esta incompleta: ",
                           conditionMessage(e))
)
if (!identical(as.integer(dim(matriz)), c(29733L, 100064L))) {
  stop("Dimensiones inesperadas: ", paste(dim(matriz), collapse = " x "))
}
if (anyDuplicated(colnames(matriz)) || anyDuplicated(rownames(matriz))) {
  stop("Hay nombres de genes o barcodes duplicados.")
}

# 01C. Leer metadatos originales y alinearlos con las columnas ---------------
m <- utils::read.csv(ruta_meta, stringsAsFactors = FALSE,
                     check.names = FALSE, na.strings = c("", "NA"))
if (nrow(m) && identical(as.character(m$NAME[1]), "TYPE")) {
  m <- m[-1, , drop = FALSE]  # Fila tecnica del archivo del Broad Portal
}

requeridas <- c("NAME", "Patient", "celltype_major", "subtype",
                "nCount_RNA", "nFeature_RNA", "Percent_mito")
if (!all(requeridas %in% names(m))) {
  stop("Faltan columnas en metadatos: ",
       paste(setdiff(requeridas, names(m)), collapse = ", "))
}
for (col in c("nCount_RNA", "nFeature_RNA", "Percent_mito")) {
  m[[col]] <- suppressWarnings(as.numeric(m[[col]]))
  if (anyNA(m[[col]])) stop("Valores de QC no numericos o faltantes en: ", col)
}
if (anyNA(m$NAME) || anyDuplicated(m$NAME) ||
    !setequal(m$NAME, colnames(matriz))) {
  stop("Los identificadores celulares de metadatos y matriz no coinciden.")
}

m <- m[match(colnames(matriz), m$NAME), , drop = FALSE]
rownames(m) <- m$NAME
stopifnot(identical(rownames(m), colnames(matriz)))

# Conservamos indicadores publicados para poder contrastarlos en el objetivo A.
names(m)[names(m) == "nCount_RNA"]   <- "nCount_RNA_original"
names(m)[names(m) == "nFeature_RNA"] <- "nFeature_RNA_original"

# 01D. Crear objeto inicial o verificar uno ya guardado -----------------------
if (file.exists(ruta_rds)) {
  cat("[A/01] Seurat inicial ya existe; se vuelve a abrir y validar.\n")
  obj <- readRDS(ruta_rds)
} else {
  cat("[A/01] Creando objeto Seurat sin filtrar...\n")
  obj <- Seurat::CreateSeuratObject(
    counts = matriz,
    meta.data = m,
    project = "GSE176078",
    min.cells = 0,
    min.features = 0
  )
  saveRDS(obj, ruta_rds)
}

if (nrow(obj) != 29733L || ncol(obj) != 100064L ||
    !identical(colnames(obj), colnames(matriz)) ||
    !identical(rownames(obj), rownames(matriz))) {
  stop("El objeto Seurat no coincide con la matriz original.")
}
if (!all(c("Patient", "celltype_major", "subtype", "Percent_mito",
           "nCount_RNA_original", "nFeature_RNA_original") %in%
         colnames(obj[[]]))) {
  stop("Faltan anotaciones esperadas en el objeto Seurat.")
}

conteos <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
if (!inherits(conteos, "IterableMatrix")) {
  # En BPCells 0.3.1 Seurat puede mostrar la clase concreta "RenameDims".
  if (!identical(attr(class(conteos), "package"), "BPCells")) {
    stop("La capa counts no parece estar respaldada en disco por BPCells.")
  }
}
# Leer una submatriz de 2 x 2 confirma el acceso real a BPCells.
# No usar as.numeric() directamente sobre IterableMatrix (objeto S4).
# Convertimos SOLO cuatro valores, nunca la matriz completa.
tryCatch({
  muestra <- as.matrix(conteos[1:2, 1:2])
  stopifnot(is.matrix(muestra), all(is.finite(muestra)))
}, error = function(e) {
  stop("No se pudo leer una submatriz BPCells: ", conditionMessage(e))
})

# 01E. Registrar resultados y versiones para reproducibilidad ----------------
resumen <- data.frame(
  estudio = "GSE176078",
  genes = nrow(obj),
  celulas = ncol(obj),
  donantes = length(unique(obj$Patient)),
  capa = paste(SeuratObject::Layers(obj[["RNA"]]), collapse = ";"),
  tipo_conteos = paste(class(conteos), collapse = ";"),
  filtros_aplicados = "ninguno",
  normalizacion = "ninguna",
  ruta_datos_procesados = if (
    nzchar(Sys.getenv("RNASEQ_DATA_DIR", unset = ""))
  ) {
    "RNASEQ_DATA_DIR/GSE176078"
  } else {
    "data/processed/GSE176078"
  }
)
utils::write.csv(resumen, "docs/A_Analisis_SingleCell/01_objeto_seurat/resumen.csv", row.names = FALSE)
versiones <- data.frame(
  componente = c("R", paquetes),
  version = c(as.character(getRversion()),
              vapply(paquetes, function(p) as.character(utils::packageVersion(p)), ""))
)
utils::write.csv(versiones, "docs/A_Analisis_SingleCell/01_objeto_seurat/versiones.csv", row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()), "logs/A_Analisis_SingleCell/01_objeto_seurat/session_info.txt")

cat("\n== OBJETIVO A / ETAPA 01 COMPLETADA ==\n")
print(resumen[, c("genes", "celulas", "donantes", "capa", "tipo_conteos")],
      row.names = FALSE)
cat("Objeto guardado en:", ruta_rds, "\n")
cat("Registro: docs/A_Analisis_SingleCell/01_objeto_seurat/ y logs/A_Analisis_SingleCell/01_objeto_seurat/session_info.txt\n")
cat("Sin filtros, sin normalizacion, sin PCA y sin UMAP.\n")
