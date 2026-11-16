# ============================================================
# BIOINFORMATICA II - Proyecto RNA-seq bulk
# Objetivo A / etapa 03, QC prefiltrado por tipo celular, donante y GENES
# Entrada: data/processed/GSE176078/seurat_qc_inicial.rds
# Salidas: docs/A_Analisis_SingleCell/03_qc_por_grupos/, figures/A_Analisis_SingleCell/03_qc_por_grupos/, logs/
# Calcula deteccion de genes con BPCells sin cargar la matriz completa a RAM.
# No aplica filtros, no normaliza y no altera el objeto Seurat.
# Ejecutar desde la raiz del proyecto:
# source("scripts/A_Analisis_SingleCell/03_comparar_qc_grupos.R")
# ============================================================

if (!file.exists("scripts/A_Analisis_SingleCell/02_qc_inicial.R")) {
  stop("Abre la carpeta raiz del proyecto y guarda primero el script A/02.")
}
for (paquete in c("ggplot2", "BPCells", "SeuratObject")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("Falta el paquete: ", paquete)
  }
}

base_datos <- Sys.getenv("RNASEQ_DATA_DIR", unset = "")
if (!nzchar(base_datos)) base_datos <- file.path("data", "processed")
entrada <- file.path(base_datos, "GSE176078", "seurat_qc_inicial.rds")
if (!file.exists(entrada)) {
  stop("No se encuentra el objeto con percent.mt: ", entrada,
       "\nEjecuta antes scripts/A_Analisis_SingleCell/02_qc_inicial.R")
}

salida_docs <- "docs/A_Analisis_SingleCell/03_qc_por_grupos"
salida_figuras <- "figures/A_Analisis_SingleCell/03_qc_por_grupos"
dir.create(salida_docs, recursive = TRUE, showWarnings = FALSE)
dir.create(salida_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create("logs/A_Analisis_SingleCell/03_qc_por_grupos", recursive = TRUE, showWarnings = FALSE)

# 1. Abrir objeto; SOLO usamos sus metadatos para la comparacion.
obj <- readRDS(entrada)
meta <- obj[[]]
stopifnot(nrow(obj) == 29733L, ncol(obj) == 100064L,
          length(unique(meta$Patient)) == 26L)
indicadores <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
agrupadores <- c("celltype_major", "Patient")
necesarias <- c(indicadores, agrupadores)
if (!all(necesarias %in% colnames(meta))) {
  stop("Faltan columnas: ", paste(setdiff(necesarias, colnames(meta)), collapse = ", "))
}
if (!identical(rownames(meta), colnames(obj))) {
  stop("El orden de los metadatos no coincide con el de las celulas.")
}
for (nm in indicadores) {
  if (!is.numeric(meta[[nm]]) || anyNA(meta[[nm]]) ||
      any(!is.finite(meta[[nm]]))) {
    stop("Metrica QC invalida: ", nm)
  }
}
for (nm in agrupadores) {
  if (anyNA(meta[[nm]]) || any(!nzchar(as.character(meta[[nm]])))) {
    stop("Anotaciones ausentes en: ", nm)
  }
}
if (any(meta$nCount_RNA <= 0)) stop("Se encontraron UMI no positivos.")

cat("\n== QC POR GRUPOS, ANTES DE FILTRAR ==\n")
cat("Celulas:", nrow(meta), "\n")
cat("Tipos celulares:", length(unique(meta$celltype_major)), "\n")
cat("Donantes:", length(unique(meta$Patient)), "\n")

# 2. Resumir los tres indicadores por tipo celular y por donante.
# Los percentiles se calculan sobre CELULAS, no sobre donantes.
resumir_por <- function(datos, columna_grupo) {
  indices <- split(seq_len(nrow(datos)), as.character(datos[[columna_grupo]]))
  tablas <- lapply(names(indices), function(nombre_grupo) {
    ix <- indices[[nombre_grupo]]
    do.call(rbind, lapply(indicadores, function(indicador) {
      x <- datos[[indicador]][ix]
      q <- unname(stats::quantile(x, probs = c(0.25, 0.5, 0.75)))
      data.frame(
        grupo = nombre_grupo,
        indicador = indicador,
        n_celulas = length(ix),
        minimo = min(x),
        p25 = q[1],
        mediana = q[2],
        p75 = q[3],
        maximo = max(x),
        stringsAsFactors = FALSE
      )
    }))
  })
  tabla <- do.call(rbind, tablas)
  rownames(tabla) <- NULL
  names(tabla)[1] <- columna_grupo
  tabla
}

por_tipo <- resumir_por(meta, "celltype_major")
por_donante <- resumir_por(meta, "Patient")
utils::write.csv(por_tipo, file.path(salida_docs, "resumen_qc_por_tipo_pre.csv"),
                 row.names = FALSE)
utils::write.csv(por_donante, file.path(salida_docs, "resumen_qc_por_donante_pre.csv"),
                 row.names = FALSE)

# 3. Contar celulas por combinacion de donante y tipo.
# Incluye combinaciones con cero celulas: relevante para referencia posterior.
n_por_donante_tipo <- as.data.frame(table(
  Patient = meta$Patient, celltype_major = meta$celltype_major
), stringsAsFactors = FALSE)
names(n_por_donante_tipo)[3] <- "n_celulas"
stopifnot(sum(n_por_donante_tipo$n_celulas) == nrow(meta))
utils::write.csv(n_por_donante_tipo,
                 file.path(salida_docs, "celulas_por_donante_y_tipo_pre.csv"),
                 row.names = FALSE)

# 4. Graficos por grupo: 3 indicadores en paneles separados.
# UMI se presenta como log10 para evitar que los extremos oculten el resto.
# En los boxplots ocultamos los PUNTOS atipicos para legibilidad,
# pero NO eliminamos esas celulas; minimos y maximos quedan en las tablas.
graficar_qc <- function(datos, columna_grupo, titulo, archivo, ancho) {
  grupos <- factor(as.character(datos[[columna_grupo]]),
                   levels = sort(unique(as.character(datos[[columna_grupo]]))))
  n <- nrow(datos)
  long <- data.frame(
    grupo = rep(grupos, times = 3),
    indicador = factor(
      rep(c("Genes detectados (nFeature_RNA)",
            "Conteos UMI (log10)",
            "Mitocondrial (%)"), each = n),
      levels = c("Genes detectados (nFeature_RNA)",
                 "Conteos UMI (log10)", "Mitocondrial (%)")
    ),
    valor = c(datos$nFeature_RNA, log10(datos$nCount_RNA), datos[["percent.mt"]])
  )
  grafico <- ggplot2::ggplot(long, ggplot2::aes(x = grupo, y = valor)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.7,
                          fill = "#A7C9DA", color = "#365C72") +
    ggplot2::facet_grid(indicador ~ ., scales = "free_y") +
    ggplot2::labs(title = titulo,
                  subtitle = "Antes del filtrado; los puntos extremos se omiten solo del dibujo",
                  x = NULL, y = "Valor del indicador") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 60, hjust = 1, vjust = 1))
  ggplot2::ggsave(file.path(salida_figuras, archivo), plot = grafico,
                  width = ancho, height = 10, dpi = 150)
  invisible(NULL)
}

graficar_qc(meta, "celltype_major", "QC por tipo celular - GSE176078",
            "qc_por_tipo_celular_pre.png", ancho = 12)
graficar_qc(meta, "Patient", "QC por donante - GSE176078",
            "qc_por_donante_pre.png", ancho = 15)

# 5. Mapa de numeros de celulas por donante y tipo.
heatmap <- ggplot2::ggplot(n_por_donante_tipo,
                          ggplot2::aes(x = celltype_major, y = Patient,
                                       fill = n_celulas)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::scale_fill_viridis_c(name = "Celulas") +
  ggplot2::labs(title = "Celulas por donante y tipo celular (antes de filtrar)",
                x = "Tipo celular", y = "Donante") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 45, hjust = 1))
ggplot2::ggsave(file.path(salida_figuras, "celulas_por_donante_tipo_pre.png"),
                plot = heatmap, width = 12, height = 10, dpi = 150)

# 6. QC DE GENES (integrado al paso A/03, SIN filtrado) ------------
# Para cada gen, medimos EN CUANTAS CELULAS hay al menos 1 UMI.
# BPCells::matrix_stats procesa la matriz por bloques en disco.
# NO usar as.matrix(matriz_counts): excederia la RAM disponible.
cat("\n== QC DE GENES - DETECCION CELULAR ==\n")

matriz_counts <- SeuratObject::LayerData(
  obj, assay = "RNA", layer = "counts"
)
stopifnot(
  nrow(matriz_counts) == nrow(obj),
  ncol(matriz_counts) == ncol(obj),
  identical(rownames(matriz_counts), rownames(obj)),
  identical(colnames(matriz_counts), colnames(obj))
)

# Una lectura secuenciada y un UMI son conceptos distintos.
# Aqui evaluamos UMI crudos, no lecturas FASTQ.
resultado_stats <- BPCells::matrix_stats(
  matriz_counts,
  row_stats = "nonzero",
  col_stats = "none",
  threads = 1L
)
estadisticas_filas <- resultado_stats$row_stats
if (!is.matrix(estadisticas_filas) ||
    ncol(estadisticas_filas) != nrow(obj)) {
  stop("BPCells devolvio estadisticas de genes con dimensiones inesperadas.")
}
indice_no_cero <- if ("nonzero" %in% rownames(estadisticas_filas)) {
  match("nonzero", rownames(estadisticas_filas))
} else {
  1L
}
n_celulas_gen <- as.numeric(estadisticas_filas[indice_no_cero, ])
if (length(n_celulas_gen) != nrow(obj) ||
    anyNA(n_celulas_gen) ||
    any(n_celulas_gen < 0 | n_celulas_gen > ncol(obj))) {
  stop("Conteo de celulas por gen invalido.")
}

qc_genes <- data.frame(
  gen = rownames(obj),
  n_celulas_detectado = n_celulas_gen,
  porcentaje_celulas = 100 * n_celulas_gen / ncol(obj),
  stringsAsFactors = FALSE
)
stopifnot(!anyDuplicated(qc_genes$gen))
utils::write.csv(
  qc_genes,
  file.path(salida_docs, "qc_genes_deteccion_por_gen_pre.csv"),
  row.names = FALSE
)

# Las cifras de deteccion se refieren a genes, no a celulas.
# Los cortes 1, 3 y 10 son ESCENARIOS; no se aplican aqui.
umbrales_genes <- c(1L, 3L, 10L)
impacto_genes <- data.frame(
  umbral_min_celulas = umbrales_genes,
  genes_conservados = vapply(
    umbrales_genes, function(u) sum(n_celulas_gen >= u), numeric(1)
  ),
  stringsAsFactors = FALSE
)
impacto_genes$genes_excluidos_hipoteticos <- nrow(obj) - impacto_genes$genes_conservados
impacto_genes$porcentaje_excluido_hipotetico <- round(
  100 * impacto_genes$genes_excluidos_hipoteticos / nrow(obj), 2
)
utils::write.csv(
  impacto_genes,
  file.path(salida_docs, "qc_genes_impacto_umbrales_pre.csv"),
  row.names = FALSE
)

resumen_genes <- data.frame(
  genes_totales = nrow(obj),
  genes_no_detectados = sum(n_celulas_gen == 0),
  genes_detectados_en_una_celula = sum(n_celulas_gen == 1),
  genes_detectados_en_dos_celulas = sum(n_celulas_gen == 2),
  mediana_celulas_por_gen = stats::median(n_celulas_gen),
  max_celulas_por_gen = max(n_celulas_gen),
  stringsAsFactors = FALSE
)
utils::write.csv(
  resumen_genes,
  file.path(salida_docs, "qc_genes_resumen_pre.csv"),
  row.names = FALSE
)

p_deteccion <- ggplot2::ggplot(
  qc_genes,
  ggplot2::aes(x = log10(n_celulas_detectado + 1))
) +
  ggplot2::geom_histogram(bins = 70, fill = "#679CB7", color = "white") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "QC de genes: deteccion en celulas - GSE176078",
    subtitle = "Conteos crudos UMI, antes de cualquier filtro adicional",
    x = "log10(numero de celulas que detectan el gen + 1)",
    y = "Numero de genes"
  )
ggplot2::ggsave(
  file.path(salida_figuras, "qc_genes_distribucion_deteccion_pre.png"),
  plot = p_deteccion, width = 9, height = 6, dpi = 160
)

impacto_largo <- rbind(
  data.frame(umbral = factor(paste0(">= ", umbrales_genes, " celulas"),
                               levels = paste0(">= ", umbrales_genes, " celulas")),
             resultado = "Conservados", n_genes = impacto_genes$genes_conservados),
  data.frame(umbral = factor(paste0(">= ", umbrales_genes, " celulas"),
                               levels = paste0(">= ", umbrales_genes, " celulas")),
             resultado = "Excluidos (hipotetico)",
             n_genes = impacto_genes$genes_excluidos_hipoteticos)
)
p_impacto <- ggplot2::ggplot(
  impacto_largo,
  ggplot2::aes(x = umbral, y = n_genes, fill = resultado)
) +
  ggplot2::geom_col() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Impacto HIPOTETICO de criterios de deteccion por gen",
    subtitle = "Solo simulacion; no se ha eliminado ningun gen",
    x = "Criterio de deteccion", y = "Numero de genes", fill = "Estado"
  )
ggplot2::ggsave(
  file.path(salida_figuras, "qc_genes_escenarios_pre.png"),
  plot = p_impacto, width = 9, height = 6, dpi = 160
)

cat("Genes totales:", nrow(qc_genes), "\n")
cat("Genes no detectados:", resumen_genes$genes_no_detectados, "\n")
cat("Mediana de celulas por gen:", resumen_genes$mediana_celulas_por_gen, "\n")
cat("\nImpacto hipotetico de criterios de genes:\n")
print(impacto_genes, row.names = FALSE)

# Registrar interpretacion y parametros sin alterar el objeto Seurat.
writeLines(c(
  "A/03 - QC de genes: deteccion global en la matriz sin filtrar.",
  "Un gen esta detectado en una celula si cuenta con al menos 1 UMI.",
  "Los umbrales 1, 3 y 10 celulas son simulaciones, NO filtros aplicados.",
  "No se ha escogido un umbral definitivo ni se han excluido genes o celulas.",
  "El QC celular por tipo y donante esta en los otros resultados del paso 03.",
  "La seleccion de marcadores para la referencia corresponde al objetivo B, tras separar donantes."
), file.path(salida_docs, "qc_genes_notas_metodologicas.txt"))

# Liberar objetos intermedios sin alterar el objeto Seurat original.
rm(resultado_stats, estadisticas_filas, matriz_counts)

# Guardar decisiones de visualizacion para reproducibilidad.
parametros <- data.frame(
  parametro = c("datos", "filtros", "indicadores", "UMI_grafico", "outliers_boxplot",
                "unidad_resumen_celular", "criterio_deteccion_gen", "umbrales_explorados_genes"),
  valor = c("seurat_qc_inicial.rds", "ninguno", paste(indicadores, collapse = ";"),
            "log10(nCount_RNA)", "ocultos solo visualmente", "celula",
            ">= 1 UMI por gen en celula", "1;3;10 celulas (solo simulacion)")
)
utils::write.csv(parametros,
                 file.path(salida_docs, "parametros_qc_grupos.csv"), row.names = FALSE)

# 6. Registrar el entorno de software.
writeLines(capture.output(utils::sessionInfo()),
           "logs/A_Analisis_SingleCell/03_qc_por_grupos/session_info.txt")
cat("\n== ARCHIVOS GENERADOS ==\n")
print(list.files(salida_docs, full.names = TRUE))
print(list.files(salida_figuras, full.names = TRUE))
cat("\nPASO A/03 - QC CELULAR POR GRUPOS Y QC DE GENES COMPLETADO.\n")
cat("NO se aplicaron filtros a celulas ni genes; no se modificaron los conteos.\n")

# ============================================================
# A/03 - QC GENES: EXPRESION Y FRECUENCIA DE DETECCION
# ============================================================
# No se eliminan genes ni celulas.
# Se utilizan conteos UMI crudos almacenados en BPCells.

cat("\n== QC GENES: EXPRESION Y DETECCION ==\n")

# 1. Obtener estadisticas por gen desde BPCells
matriz_expr <- SeuratObject::LayerData(
  obj,
  assay = "RNA",
  layer = "counts"
)

estadisticas <- BPCells::matrix_stats(
  matriz_expr,
  row_stats = "mean",
  col_stats = "none",
  threads = 1L
)$row_stats

stopifnot(
  all(c("nonzero", "mean") %in% rownames(estadisticas)),
  ncol(estadisticas) == nrow(obj)
)

# 2. Preparar tabla de QC por gen
qc_genes_extra <- data.frame(
  gen = rownames(obj),
  n_celulas = as.numeric(estadisticas["nonzero", ]),
  umi_medio = as.numeric(estadisticas["mean", ])
)

qc_genes_extra$porcentaje_deteccion <-
  100 * qc_genes_extra$n_celulas / ncol(obj)

qc_genes_extra$umi_total <- round(
  qc_genes_extra$umi_medio * ncol(obj)
)

# 3. Clasificar genes para los graficos
qc_genes_extra$categoria <- ifelse(
  grepl("^MT-", qc_genes_extra$gen),
  "Mitocondrial",
  ifelse(
    grepl("^RP[LS]", qc_genes_extra$gen),
    "Ribosomal",
    "Otros"
  )
)

# Validar contra UMI totales del objeto Seurat
diferencia_relativa <- abs(
  sum(qc_genes_extra$umi_total) - sum(obj$nCount_RNA)
) / sum(obj$nCount_RNA)

stopifnot(diferencia_relativa < 1e-5)

# Guardar tabla completa
write.csv(
  qc_genes_extra,
  file.path(salida_docs, "qc_genes_exp_media_y_deteccion_pre.csv"),
  row.names = FALSE
)

# ============================================================
# GRAFICO 1: EXPRESION MEDIA VS DETECCION
# ============================================================

p_expresion <- ggplot2::ggplot(
  qc_genes_extra,
  ggplot2::aes(
    x = log10(umi_medio + 0.0001),
    y = porcentaje_deteccion,
    color = categoria
  )
) +
  ggplot2::geom_point(
    alpha = 0.4,
    size = 0.8
  ) +
  ggplot2::labs(
    title = "QC genes - Expresion media vs deteccion",
    subtitle = "GSE176078 | Antes de filtrar",
    x = "log10(UMI promedio por celula + 0.0001)",
    y = "Celulas que detectan el gen (%)",
    color = "Categoria"
  ) +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(
  file.path(
    salida_figuras,
    "qc_genes_expresion_vs_deteccion_pre.png"
  ),
  plot = p_expresion,
  width = 10,
  height = 6,
  dpi = 160,
  bg = "white"
)

# ============================================================
# GRAFICO 2: TOP 20 GENES CON MAS UMI
# ============================================================

top20 <- qc_genes_extra[
  order(-qc_genes_extra$umi_total, qc_genes_extra$gen),
]

top20 <- head(top20, 20)

write.csv(
  top20,
  file.path(salida_docs, "qc_genes_top20_umi_pre.csv"),
  row.names = FALSE
)

p_top20 <- ggplot2::ggplot(
  top20,
  ggplot2::aes(
    x = reorder(gen, umi_total),
    y = umi_total / 1000000,
    fill = categoria
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "QC genes - Top 20 por conteos UMI",
    subtitle = "GSE176078 | Conteos crudos acumulados",
    x = "Gen",
    y = "Millones de UMI",
    fill = "Categoria"
  ) +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(
  file.path(
    salida_figuras,
    "qc_genes_top20_umi_pre.png"
  ),
  plot = p_top20,
  width = 10,
  height = 7,
  dpi = 160,
  bg = "white"
)

# Liberar objetos temporales
rm(matriz_expr, estadisticas)

cat("\n== RESULTADOS QC DE GENES ==\n")
cat("Genes evaluados:", nrow(qc_genes_extra), "\n")
cat("Diferencia relativa de UMI:", diferencia_relativa, "\n")
cat("Graficos adicionales guardados correctamente.\n")
cat("No se aplicaron filtros.\n")