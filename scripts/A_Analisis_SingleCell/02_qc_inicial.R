# =============================================================
# BIOINFORMATICA II - PROYECTO RNA-seq bulk
# OBJETIVO A / ETAPA 02 - QC inicial de celulas (ANTES de filtrar)
# Entrada: data/processed/GSE176078/seurat_inicial.rds
# Salidas: docs/A_Analisis_SingleCell/02_qc_inicial/ y figures/A_Analisis_SingleCell/02_qc_inicial/
# Calcula percent.mt desde los conteos crudos con BPCells.
# Este script NO filtra, normaliza ni modifica el objeto original.
# Ejecutar desde la raiz del proyecto:
# source("scripts/A_Analisis_SingleCell/02_qc_inicial.R")
# =============================================================

if (!file.exists("scripts/A_Analisis_SingleCell/01_crear_objeto_seurat.R")) {
  stop("Abre la carpeta raiz 'Proyecto RNA-seq bulk' en Positron.")
}

for (paquete in c("SeuratObject", "ggplot2")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("Falta instalar el paquete: ", paquete)
  }
}

# Consistente con la ubicación opcional del paso 02.
base_datos <- Sys.getenv("RNASEQ_DATA_DIR", unset = "")
if (!nzchar(base_datos)) base_datos <- file.path("data", "processed")
carpeta_objeto <- file.path(base_datos, "GSE176078")
entrada <- file.path(carpeta_objeto, "seurat_inicial.rds")
if (!file.exists(entrada)) stop("No encuentro el objeto: ", entrada)

dir.create("docs/A_Analisis_SingleCell/02_qc_inicial", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/A_Analisis_SingleCell/02_qc_inicial", recursive = TRUE, showWarnings = FALSE)
dir.create("logs/A_Analisis_SingleCell/02_qc_inicial", recursive = TRUE, showWarnings = FALSE)

# 1. Abrir objeto Seurat. Los conteos permanecen en BPCells en disco.
obj <- readRDS(entrada)
stopifnot(nrow(obj) == 29733L, ncol(obj) == 100064L)
meta <- obj[[]]
columnas <- c("nFeature_RNA", "nCount_RNA", "Percent_mito",
              "Patient", "celltype_major")
if (!all(columnas %in% colnames(meta))) {
  stop("Faltan columnas requeridas: ",
       paste(setdiff(columnas, colnames(meta)), collapse = ", "))
}
for (nombre in columnas[1:3]) {
  if (!is.numeric(meta[[nombre]]) || anyNA(meta[[nombre]])) {
    stop("Indicador de QC no numerico o con NA: ", nombre)
  }
}

# 2. Identificar genes mitocondriales humanos y calcular el porcentaje
#    desde los conteos UMI crudos, sin cargar la matriz entera en RAM.
#    Seurat usa el prefijo ^MT- para genes mitocondriales humanos.
#    Convertimos SOLO las filas MT a matriz convencional (~10 MB),
#    nunca toda la matriz; esto evita problemas de despacho colSums().
mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
if (!length(mt_genes)) {
  stop("No se encontraron genes con prefijo MT-. Revisar anotacion.")
}
utils::write.csv(
  data.frame(gen = mt_genes),
  "docs/A_Analisis_SingleCell/02_qc_inicial/genes_mitocondriales.csv", row.names = FALSE
)

matriz_counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
stopifnot(identical(rownames(meta), colnames(obj)))
# Extraer SOLO genes mitocondriales (~13 x 100064 valores): RAM controlada
conteos_mt <- as.matrix(matriz_counts[mt_genes, , drop = FALSE])

# Validacion informativa: comparar VALORES de dimensiones, no sus atributos
# ni el tipo entero/double (identical puede producir falsos negativos).
dim_esperada <- c(length(mt_genes), ncol(obj))
dim_obtenida <- dim(conteos_mt)
cat("\n== DIMENSIONES SUBMATRIZ MITOCONDRIAL ==\n")
cat("Obtenida:", paste(dim_obtenida, collapse = " x "), "\n")
cat("Esperada:", paste(dim_esperada, collapse = " x "), "\n")
if (length(dim_obtenida) != 2L ||
    !all(as.numeric(dim_obtenida) == as.numeric(dim_esperada))) {
  stop("La submatriz mitocondrial no tiene las dimensiones esperadas. ",
       "No se calculara percent.mt hasta revisar esta diferencia.")
}

# Las columnas deben seguir el orden de celulas del objeto Seurat.
# Si la conversion no conserva los nombres, el indexado original sí
# mantiene el orden; verificamos los nombres cuando vienen presentes.
if (!is.null(colnames(conteos_mt)) &&
    !identical(colnames(conteos_mt), colnames(obj))) {
  stop("El orden o los nombres celulares no coinciden con el objeto Seurat.")
}
if (!is.null(rownames(conteos_mt)) &&
    !identical(rownames(conteos_mt), mt_genes)) {
  stop("El orden o los nombres de genes MT no coinciden.")
}

umi_mt <- unname(base::colSums(conteos_mt))
rm(conteos_mt)
stopifnot(length(umi_mt) == ncol(obj), all(meta$nCount_RNA > 0))
percent_mt <- 100 * as.numeric(umi_mt) / as.numeric(meta$nCount_RNA)
stopifnot(all(is.finite(percent_mt)))

# No sobrescribir el dato publicado; guardar el recalculado con otro nombre.
obj[["percent.mt"]] <- percent_mt
meta <- obj[[]]

# Comparar ambos porcentajes celula por celula.
comparacion <- data.frame(
  NAME = colnames(obj),
  Percent_mito_publicado = meta$Percent_mito,
  percent_mt_calculado = meta[["percent.mt"]],
  diferencia = meta[["percent.mt"]] - meta$Percent_mito
)
utils::write.csv(
  comparacion, "docs/A_Analisis_SingleCell/02_qc_inicial/comparacion_mitocondrial.csv",
  row.names = FALSE
)
resumen_mt <- data.frame(
  n_genes_mt = length(mt_genes),
  mediana_mt_calculada = median(meta[["percent.mt"]]),
  mediana_mt_publicada = median(meta$Percent_mito),
  maxima_diferencia_absoluta = max(abs(comparacion$diferencia)),
  celulas_diferencia_mayor_0_01 = sum(abs(comparacion$diferencia) > 0.01)
)
utils::write.csv(
  resumen_mt, "docs/A_Analisis_SingleCell/02_qc_inicial/resumen_mitocondrial.csv",
  row.names = FALSE
)

# Objeto enriquecido para el proximo paso (el original no se modifica).
salida_obj <- file.path(carpeta_objeto, "seurat_qc_inicial.rds")
saveRDS(obj, salida_obj)

# 3. Resumen del QC sin filtrar.
resumir <- function(x, indicador) {
  data.frame(indicador = indicador, minimo = min(x),
             p25 = unname(quantile(x, 0.25)), mediana = median(x),
             media = mean(x), p75 = unname(quantile(x, 0.75)),
             maximo = max(x))
}
indicadores_qc <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
resumen <- do.call(rbind, lapply(indicadores_qc, function(nm) {
  resumir(meta[[nm]], nm)
}))
utils::write.csv(resumen, "docs/A_Analisis_SingleCell/02_qc_inicial/resumen_qc_pre.csv",
                 row.names = FALSE)

# 4. Numero de celulas por donante y tipo celular: antes de filtrar.
conteos <- as.data.frame(table(Patient = meta$Patient,
                               celltype_major = meta$celltype_major),
                         stringsAsFactors = FALSE)
names(conteos)[3] <- "n_celulas"
utils::write.csv(conteos,
                 "docs/A_Analisis_SingleCell/02_qc_inicial/celulas_por_donante_y_tipo_pre.csv",
                 row.names = FALSE)

# 5. Violines previos al filtrado con el percent.mt recien calculado.
violines <- rbind(
  data.frame(indicador = "Genes detectados (nFeature_RNA)",
             valor = meta$nFeature_RNA),
  data.frame(indicador = "Conteos UMI (nCount_RNA)",
             valor = meta$nCount_RNA),
  data.frame(indicador = "Mitocondrial (%) [calculado]",
             valor = meta[["percent.mt"]])
)
violines$indicador <- factor(violines$indicador,
  levels = c("Genes detectados (nFeature_RNA)",
             "Conteos UMI (nCount_RNA)",
             "Mitocondrial (%) [calculado]"))

p_violin <- ggplot2::ggplot(violines,
                            ggplot2::aes(x = indicador, y = valor)) +
  ggplot2::geom_violin(fill = "#97BDD2", color = "#43667A", trim = TRUE) +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA,
                        fill = "white", alpha = 0.9) +
  ggplot2::facet_wrap(~ indicador, scales = "free", nrow = 1) +
  ggplot2::labs(title = "QC inicial - antes de aplicar filtros",
                x = NULL, y = "Valor del indicador") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())

ggplot2::ggsave("figures/A_Analisis_SingleCell/02_qc_inicial/violines_qc_pre.png",
                plot = p_violin, width = 12, height = 5, dpi = 160)

# 6. Relacion entre genes detectados y conteos UMI.
p_relacion <- ggplot2::ggplot(meta,
   ggplot2::aes(x = nCount_RNA, y = nFeature_RNA)) +
  ggplot2::geom_point(size = 0.25, alpha = 0.13, color = "#336B86") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(title = "QC inicial - genes detectados vs conteos UMI",
                x = "Conteos UMI (escala log10)", y = "Genes detectados") +
  ggplot2::theme_bw(base_size = 11)

ggplot2::ggsave("figures/A_Analisis_SingleCell/02_qc_inicial/genes_vs_umi_pre.png",
                plot = p_relacion, width = 7, height = 5, dpi = 160)


# 6B. Grafico FeatureScatter de Seurat, como en la clase de Noemi.
# Usa metadatos QC; no extrae ni convierte la matriz completa BPCells.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Falta instalar el paquete Seurat para FeatureScatter.")
}
p_scatter <- Seurat::FeatureScatter(
  object = obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 0.1,
  log = TRUE,
  plot.cor = FALSE,
  raster = FALSE
) +
  # Seurat asigna por defecto colores a las identidades activas (26 donantes).
  # En este grafico de QC interesa la relacion general; usamos un solo tono.
  ggplot2::scale_color_manual(
    values = rep("#35769A", nlevels(SeuratObject::Idents(obj)))
  ) +
  ggplot2::labs(
    title = "QC inicial - FeatureScatter de Seurat",
    x = "Conteos UMI (escala log10)",
    y = "Genes detectados (nFeature_RNA)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )
ggplot2::ggsave(
  "figures/A_Analisis_SingleCell/02_qc_inicial/featurescatter_pre.png",
  plot = p_scatter,
  width = 9, height = 6, dpi = 160, bg = "white"
)

# 7. Registro del entorno para reproducibilidad.
writeLines(capture.output(utils::sessionInfo()),
           "logs/A_Analisis_SingleCell/02_qc_inicial/session_info.txt")

cat("\n== OBJETIVO A1: QC INICIAL SIN FILTROS ==\n")
cat("Genes:", nrow(obj), "\n")
cat("Celulas:", ncol(obj), "\n")
cat("Donantes:", length(unique(meta$Patient)), "\n")
cat("Genes mitocondriales identificados:", length(mt_genes), "\n")
print(mt_genes)
cat("\n== COMPARACION DE PORCENTAJE MITOCONDRIAL ==\n")
print(resumen_mt, row.names = FALSE)
print(resumen, row.names = FALSE)
cat("Objeto QC guardado en:", salida_obj, "\n")
cat("\nArchivos guardados en docs/A_Analisis_SingleCell/02_qc_inicial y figures/A_Analisis_SingleCell/02_qc_inicial.\n")
cat("NO se filtraron celulas ni genes; NO se normalizo.\n")
cat("PASO A1 - GRAFICOS INICIALES COMPLETADO.\n")
