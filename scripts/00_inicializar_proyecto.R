# =============================================================
# BIOINFORMATICA II - PROYECTO RNA-seq bulk
# 00 - Inicializacion reproducible y estructura por objetivos
# =============================================================
# Ejecutar desde la raiz del repositorio:
# source("scripts/00_inicializar_proyecto.R")
#
# Este script se puede ejecutar varias veces.
# No elimina, mueve ni sobrescribe datos o resultados existentes.

# 1. Verificar la carpeta de trabajo ----------------------------
raiz <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (basename(raiz) != "Proyecto RNA-seq bulk") {
  stop("Abre primero la carpeta 'Proyecto RNA-seq bulk' en Positron.")
}

# 2. Definir los seis objetivos del trabajo ---------------------
objetivos <- c(
  "A_Analisis_SingleCell",
  "B_Referencia_Celular",
  "C_Pseudobulk_Deconvolucion",
  "D_Evaluacion_Modelo",
  "E_Analisis_Bulk",
  "F_Discusion_Critica"
)

# Etapas del objetivo A (los scripts .R iran en scripts/A_.../)
etapas_A <- c(
  "01_objeto_seurat",
  "02_qc_inicial",
  "03_qc_por_grupos",
  "04_filtrado_qc",
  "05_normalizacion_hvg",
  "06_pca",
  "07_clustering_umap",
  "08_anotacion_celular"
)

# 3. Crear carpetas con rutas relativas -------------------------
carpetas_base <- c(
  "data/raw/GSE176078",
  "data/raw/TCGA_BRCA",
  "data/processed/GSE176078",
  "scripts", "docs", "figures", "logs"
)

carpetas_objetivos <- unlist(
  lapply(c("scripts", "docs", "figures"), function(carpeta) {
    file.path(carpeta, objetivos)
  }),
  use.names = FALSE
)

# En docs y figures, las etapas A tendran su propia subcarpeta.
carpetas_etapas_A <- unlist(
  lapply(c("docs", "figures"), function(carpeta) {
    file.path(carpeta, objetivos[1], etapas_A)
  }),
  use.names = FALSE
)

carpetas <- unique(c(
  carpetas_base,
  carpetas_objetivos,
  carpetas_etapas_A
))

nuevas <- 0L
for (carpeta in carpetas) {
  if (!dir.exists(carpeta)) {
    dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
    nuevas <- nuevas + 1L
  }
  if (!dir.exists(carpeta)) stop("No se pudo crear: ", carpeta)
}

# Git no guarda directorios vacios; crear .gitkeep cuando corresponda.
# Los datos originales y procesados se excluiran con .gitignore.
carpetas_para_git <- c(
  carpetas_objetivos,
  carpetas_etapas_A
)
for (carpeta in carpetas_para_git) {
  contenido <- list.files(carpeta, all.files = TRUE, no.. = TRUE)
  if (length(contenido) == 0L) {
    file.create(file.path(carpeta, ".gitkeep"))
  }
}

# 4. Mantener actualizado .gitignore sin borrar su contenido ----
ignorar <- c(
  "data/raw/",
  "data/processed/",
  ".Rhistory",
  ".RData",
  ".Rproj.user/",
  "logs/",
  ".DS_Store"
)

if (file.exists(".gitignore")) {
  actual <- readLines(".gitignore", warn = FALSE, encoding = "UTF-8")
} else {
  actual <- character(0)
}

faltantes <- setdiff(ignorar, trimws(actual))
if (length(faltantes) > 0L) {
  writeLines(c(actual, faltantes), ".gitignore", useBytes = TRUE)
}

# 5. Documentar objetivos en README, sin borrar texto anterior --
encabezado <- "## Estructura de objetivos (A-F)"
texto_objetivos <- c(
  encabezado,
  "",
  "- **A_Analisis_SingleCell:** QC, normalizacion, PCA, UMAP y anotaciones.",
  "- **B_Referencia_Celular:** separacion de donantes y referencia.",
  "- **C_Pseudobulk_Deconvolucion:** mezclas y deconvolucion.",
  "- **D_Evaluacion_Modelo:** MAE, RMSE y comparaciones.",
  "- **E_Analisis_Bulk:** aplicacion a muestras bulk reales.",
  "- **F_Discusion_Critica:** limitaciones y analisis critico.",
  "",
  "Los scripts estan en scripts/<objetivo>/; las tablas en docs/<objetivo>/",
  "y las figuras en figures/<objetivo>/. Los datos grandes no se suben a GitHub."
)

if (!file.exists("README.md")) {
  writeLines(c(
    "# Proyecto RNA-seq bulk", "",
    "Bioinformatica II - Magister en Informatica Medica", "",
    "## Datos", "",
    "- GSE176078: scRNA-seq de cancer de mama.",
    "- TCGA-BRCA: RNA-seq bulk.", "",
    texto_objetivos
  ), "README.md", useBytes = TRUE)
} else {
  readme <- readLines("README.md", warn = FALSE, encoding = "UTF-8")
  if (!(encabezado %in% readme)) {
    writeLines(c(readme, "", texto_objetivos), "README.md", useBytes = TRUE)
  }
}

# 6. Informar el resultado y alertar sobre carpetas antiguas ----
cat("\n== PROYECTO INICIALIZADO CORRECTAMENTE ==\n")
cat("Ruta:", raiz, "\n")
cat("Carpetas verificadas:", length(carpetas), "\n")
cat("Carpetas nuevas:", nuevas, "\n")
cat("Objetivos:\n")
for (objetivo in objetivos) cat(" -", objetivo, "\n")

# Por seguridad NO se renombran las carpetas A, B, C... antiguas.
antiguas <- unlist(lapply(c("scripts", "docs", "figures"), function(p) {
  file.path(p, LETTERS[1:6])
}), use.names = FALSE)
antiguas <- antiguas[dir.exists(antiguas)]
if (length(antiguas) > 0L) {
  cat("\nATENCION: hay carpetas antiguas con nombres A-F.\n")
  cat("No se movieron ni borraron; deben migrarse con cuidado.\n")
  cat(paste0(" - ", antiguas, collapse = "\n"), "\n")
}

cat("\nDatos existentes: sin modificaciones.\n")
cat("Scripts existentes: sin modificaciones.\n")
cat("PASO 00 COMPLETADO.\n")
