# =============================================================================
# Bioinformatica II | E/05 | TCGA-BRCA - Lectura y QC inicial de RNA-seq bulk
# Entradas: 50 archivos STAR Counts descargados y validados en E/04,
#           manifiesto E/03 y referencia 10x GSE176078 (solo features.tsv.gz).
# Salidas publicas: docs/E_Analisis_Bulk/05_expresion/ y figuras.
# Salidas locales: data/processed/TCGA_BRCA/05_expresion/ (excluidas de Git).
# Requiere solo R base. NO descarga, NO filtra muestras, NO deconvoluciona.
# Ejecutar en la consola R de Positron DESDE LA RAIZ del repositorio:
# source("scripts/E_Analisis_Bulk/05_leer_QC_expresion_TCGA_BRCA.R")
# =============================================================================

options(stringsAsFactors = FALSE)

raiz <- getwd()
if (!file.exists(file.path(raiz, "README.md")) ||
    !dir.exists(file.path(raiz, "scripts", "E_Analisis_Bulk"))) {
  stop("Abre la raiz de 'Proyecto RNA-seq bulk' en Positron antes de ejecutar E/05.")
}

r_manifest <- file.path("docs", "E_Analisis_Bulk", "03_piloto",
                        "03_manifest_piloto_50_FINAL.txt")
r_integridad <- file.path("docs", "E_Analisis_Bulk", "04_descarga",
                          "01_integridad_archivos.csv")
r_resumen_e04 <- file.path("docs", "E_Analisis_Bulk", "04_descarga",
                            "02_resumen_descarga.csv")
base_expresion <- file.path("data", "raw", "TCGA_BRCA", "expresion_piloto")
r_zip_ref <- file.path("data", "raw", "GSE176078", "Allcells_raw_count_out.zip")
salida_publica <- file.path("docs", "E_Analisis_Bulk", "05_expresion")
salida_fig <- file.path("figures", "E_Analisis_Bulk", "05_expresion")
salida_local <- file.path("data", "processed", "TCGA_BRCA", "05_expresion")

for (p in c(r_manifest, r_integridad, r_resumen_e04)) {
  if (!file.exists(p)) stop("Falta entrada: ", p, ". Completa antes E/03 y E/04.")
}

manifest <- read.delim(r_manifest, sep = "\t", check.names = FALSE,
                       quote = "", stringsAsFactors = FALSE)
registro <- read.csv(r_integridad, check.names = FALSE,
                     stringsAsFactors = FALSE)
res_e04 <- read.csv(r_resumen_e04, check.names = FALSE,
                    stringsAsFactors = FALSE)
requeridas_m <- c("id", "filename", "md5", "size")
requeridas_e04 <- c("file_id", "nombre", "integro", "md5_local", "size_local")
if (!all(requeridas_m %in% names(manifest)) ||
    !all(requeridas_e04 %in% names(registro))) {
  stop("Cambió el formato del manifiesto o del reporte E/04.")
}
if (nrow(manifest) != 50L || anyDuplicated(manifest$id) ||
    anyDuplicated(manifest$filename) || nrow(registro) != 50L ||
    anyDuplicated(registro$file_id)) {
  stop("Se esperaban exactamente 50 archivos unicos del piloto E/03-E/04.")
}

ix <- match(manifest$id, registro$file_id)
if (anyNA(ix)) stop("El registro E/04 no contiene todos los File ID del piloto.")
registro <- registro[ix, , drop = FALSE]
if (!identical(as.character(manifest$filename), as.character(registro$nombre)) ||
    !all(tolower(as.character(registro$integro)) == "true") ||
    !all(as.numeric(registro$size_local) == as.numeric(manifest$size)) ||
    !all(tolower(as.character(registro$md5_local)) ==
         tolower(as.character(manifest$md5)))) {
  stop("La integridad E/04 no coincide con el manifiesto final. Ejecuta E/04 de nuevo.")
}
estado_e04 <- res_e04$valor[match("Estado de descarga completa", res_e04$indicador)]
if (length(estado_e04) != 1L || is.na(estado_e04) || estado_e04 != "Verificada") {
  stop("E/04 no declara la descarga como Verificada.")
}

rutas <- file.path(base_expresion, manifest$id, manifest$filename)
if (!all(file.exists(rutas))) {
  stop("Faltan matrices en data/raw/TCGA_BRCA/expresion_piloto/. Ejecuta E/04.")
}
if (!all(as.numeric(file.info(rutas)$size) == as.numeric(manifest$size))) {
  stop("Hay tamanos de archivo distintos del manifiesto; repite E/04.")
}
# E/04 ya calculo MD5, pero comprobamos NUEVAMENTE los bytes actuales de disco.
cat("[E/05] Comprobando MD5 de las 50 matrices actuales...\n")
md5_actual <- unname(tools::md5sum(rutas))
if (anyNA(md5_actual) ||
    !all(tolower(md5_actual) == tolower(as.character(manifest$md5)))) {
  stop("Al menos un archivo cambió desde E/04. No continuar sin integridad.")
}

# Parsear solo columnas gene_id, gene_name, gene_type y unstranded. Las escalas
# TPM/FPKM y las columnas de otros modos stranded NO se usan para los conteos.
campos_esperados <- c("gene_id", "gene_name", "gene_type", "unstranded",
                     "stranded_first", "stranded_second", "tpm_unstranded",
                     "fpkm_unstranded", "fpkm_uq_unstranded")

leer_star <- function(ruta) {
  primeras <- readLines(ruta, n = 6L, warn = FALSE)
  if (!length(primeras)) stop("Archivo vacio: ", basename(ruta))
  linea_modelo <- primeras[grepl("^# gene-model:", primeras)]
  linea_header <- primeras[!grepl("^#", primeras) & nzchar(primeras)]
  if (length(linea_header) < 1L) stop("Sin encabezado STAR Counts: ", basename(ruta))
  campos <- strsplit(linea_header[[1L]], "\t", fixed = TRUE)[[1L]]
  if (!identical(campos, campos_esperados)) {
    stop("Encabezado STAR Counts inesperado en ", basename(ruta),
         ". Columnas observadas: ", paste(campos, collapse = ", "))
  }
  d <- read.delim(ruta, header = TRUE, sep = "\t", comment.char = "#",
                  quote = "", check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = c("character", "character", "character", "numeric",
                                 rep("NULL", length(campos_esperados) - 4L)))
  if (!identical(names(d), campos_esperados[1:4])) {
    stop("Error al seleccionar columnas en: ", basename(ruta))
  }
  if (anyNA(d$gene_id) || anyDuplicated(d$gene_id) ||
      anyNA(d$unstranded) || any(!is.finite(d$unstranded)) ||
      any(d$unstranded < 0) || any(d$unstranded != floor(d$unstranded))) {
    stop("IDs de genes repetidos o conteos crudos no validos en: ", basename(ruta))
  }
  list(datos = d,
       gene_model = if (length(linea_modelo)) linea_modelo[[1L]] else "No informado")
}

cat("[E/05] Leyendo 50 matrices STAR Counts: gene_id + unstranded...\n")
primer <- leer_star(rutas[[1L]])
datos_1 <- primer$datos
ids_todos <- datos_1$gene_id
tecnicas <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
if (!all(tecnicas %in% ids_todos)) {
  stop("No se encontraron las cuatro filas tecnicas esperadas de STAR.")
}
es_gen <- !(ids_todos %in% tecnicas)
if (any(grepl("^N_", ids_todos[es_gen]))) {
  stop("Hay otras filas tecnicas N_ no previstas; revisar antes de continuar.")
}
if (sum(es_gen) < 10000L) stop("Menos de 10 mil filas de genes: revisar archivos.")

genes <- ids_todos[es_gen]
anotacion <- datos_1[es_gen, c("gene_id", "gene_name", "gene_type"), drop = FALSE]
if (anyNA(anotacion$gene_name) || anyNA(anotacion$gene_type)) {
  stop("Hay anotaciones genicas faltantes; no se realizara correspondencia.")
}
modelo <- primer$gene_model
matriz <- matrix(NA_real_, nrow = length(genes), ncol = 50L,
                 dimnames = list(genes, manifest$id))
filas_metricas <- vector("list", nrow(manifest))

for (j in seq_len(nrow(manifest))) {
  if (j == 1L) {
    z <- primer
  } else {
    z <- leer_star(rutas[[j]])
  }
  d <- z$datos
  if (!identical(as.character(d$gene_id), ids_todos) ||
      !identical(as.character(d$gene_name), as.character(datos_1$gene_name)) ||
      !identical(as.character(d$gene_type), as.character(datos_1$gene_type)) ||
      !identical(z$gene_model, modelo)) {
    stop("Genes, orden, anotacion o modelo GENCODE no homogeneos en el archivo ", j,
         ". No se combinan matrices incompatibles de forma silenciosa.")
  }
  cuenta <- d$unstranded[es_gen]
  matriz[, j] <- cuenta
  categorias_tecnicas <- d$unstranded[match(tecnicas, d$gene_id)]
  filas_metricas[[j]] <- data.frame(
    file_id = manifest$id[j],
    # 'total_conteos_genicos' NO equivale a total de lecturas secuenciadas.
    total_conteos_genicos = sum(cuenta),
    genes_con_conteos_positivos = sum(cuenta > 0),
    genes_con_10_o_mas_conteos = sum(cuenta >= 10),
    porcentaje_conteos_protein_coding =
      if (sum(cuenta) > 0) 100 * sum(cuenta[anotacion$gene_type == "protein_coding"]) /
        sum(cuenta) else NA_real_,
    N_unmapped = categorias_tecnicas[1L],
    N_multimapping = categorias_tecnicas[2L],
    N_noFeature = categorias_tecnicas[3L],
    N_ambiguous = categorias_tecnicas[4L],
    stringsAsFactors = FALSE
  )
  if (j %% 10L == 0L || j == 1L) {
    cat(sprintf("  [%02d/50] Leida y validada.\n", j))
  }
  rm(d, z)
}
rm(datos_1, primer)

metricas <- do.call(rbind, filas_metricas)
if (anyNA(matriz) || any(!is.finite(matriz)) || any(matriz < 0)) {
  stop("La matriz resultante contiene valores incorrectos.")
}
if (any(metricas$total_conteos_genicos == 0)) {
  stop("Se encontro una muestra sin conteos genicos; revisar antes de seguir.")
}

# Leer EXCLUSIVAMENTE features.tsv.gz, no extraer los 534 MB del ZIP.
# En esta referencia, features.tsv.gz identifica genes principalmente por SYMBOL.
leer_referencia <- function(r_zip) {
  if (!file.exists(r_zip)) {
    return(list(estado = "Ausente: ZIP GSE176078", simbolos = character(0),
                ids_ensg = character(0), n_filas = NA_integer_, error = ""))
  }
  tmp <- tempfile("e05_ref_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  tryCatch({
    inv <- utils::unzip(r_zip, list = TRUE)$Name
    candidatas <- inv[grepl("(^|/)(features|genes)\\.tsv(\\.gz)?$", inv)]
    candidatas <- candidatas[!grepl("__MACOSX", candidatas, fixed = TRUE)]
    if (length(candidatas) != 1L) {
      stop("No hay una unica features.tsv(.gz) o genes.tsv(.gz) en el ZIP.")
    }
    utils::unzip(r_zip, files = candidatas, exdir = tmp, overwrite = TRUE)
    archivo <- file.path(tmp, candidatas[[1L]])
    if (!file.exists(archivo)) stop("No se pudo extraer features.tsv.gz.")
    f <- if (grepl("\\.gz$", archivo)) {
      read.delim(gzfile(archivo, open = "rt"), header = FALSE,
                 sep = "\t", quote = "", stringsAsFactors = FALSE)
    } else {
      read.delim(archivo, header = FALSE, sep = "\t", quote = "",
                 stringsAsFactors = FALSE)
    }
    if (ncol(f) < 2L || !nrow(f)) stop("Formato de features 10x inesperado.")
    primera <- as.character(f[[1L]])
    segunda <- as.character(f[[2L]])
    # En otras matrices 10x V1 lleva IDs ENSG y V2 simbolos. Para GSE176078
    # se documenta exactamente lo que aparezca en el archivo, sin inventar IDs.
    ids_ensg <- if (mean(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", primera)) > 0.8) {
      primera
    } else {
      character(0)
    }
    list(estado = "Leida desde ZIP GSE176078", simbolos = segunda,
         ids_ensg = ids_ensg, n_filas = nrow(f), error = "")
  }, error = function(e) {
    list(estado = "Referencia no legible", simbolos = character(0),
         ids_ensg = character(0), n_filas = NA_integer_,
         error = conditionMessage(e))
  })
}

cat("[E/05] Inspeccionando solo los nombres de genes del ZIP GSE176078...\n")
ref <- leer_referencia(r_zip_ref)

# Comparacion preliminar POR SIMBOLOS, nunca sustitucion automatica de ENSG.
simbolos_bulk <- as.character(anotacion$gene_name)
simbolos_ref <- as.character(ref$simbolos)
simbolos_bulk_validos <- simbolos_bulk[nzchar(simbolos_bulk) &
  !duplicated(simbolos_bulk) & !duplicated(simbolos_bulk, fromLast = TRUE)]
simbolos_ref_validos <- simbolos_ref[nzchar(simbolos_ref) &
  !duplicated(simbolos_ref) & !duplicated(simbolos_ref, fromLast = TRUE)]
comunes_simbolo <- intersect(simbolos_bulk_validos, simbolos_ref_validos)
# No combinar por simbolo si alguno se repite; la superposicion no es aun
# una matriz armonizada ni una estimacion de proporciones celulares.
ensg_sin_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))
id_compartidos <- if (length(ref$ids_ensg)) {
  length(intersect(ensg_sin_version(anotacion$gene_id),
                   ensg_sin_version(ref$ids_ensg)))
} else NA_integer_

# Las salidas por muestra y los objetos gene x muestra quedan solamente LOCALES.
dir.create(salida_publica, recursive = TRUE, showWarnings = FALSE)
dir.create(salida_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(salida_local, recursive = TRUE, showWarnings = FALSE)
saveRDS(matriz, file.path(salida_local, "matriz_counts_unstranded_50.rds"),
        compress = "xz")
write.csv(anotacion, file.path(salida_local, "anotacion_genes_GDC.csv"),
          row.names = FALSE)
write.csv(metricas, file.path(salida_local, "metricas_por_muestra.csv"),
          row.names = FALSE)

valores <- c(
  "Archivos STAR Counts leidos" = ncol(matriz),
  "Genes anotados en matrices bulk" = nrow(matriz),
  "Genes protein_coding anotados" = sum(anotacion$gene_type == "protein_coding"),
  "Genes con conteos positivos en alguna muestra" = sum(rowSums(matriz > 0) > 0),
  "Genes con conteos positivos en todas las muestras" = sum(rowSums(matriz > 0) == ncol(matriz)),
  "Modelos genicos distintos" = length(unique(modelo)),
  "Simbolos bulk unicos no ambiguos" = length(simbolos_bulk_validos),
  "Genes referencia GSE176078 (filas features.tsv)" = ref$n_filas,
  "Simbolos referencia unicos no ambiguos" = length(simbolos_ref_validos),
  "Simbolos compartidos no ambiguos" = length(comunes_simbolo),
  "IDs Ensembl compartidos sin version" = id_compartidos
)
resumen <- data.frame(indicador = names(valores), valor = as.numeric(valores),
                      stringsAsFactors = FALSE)
write.csv(resumen, file.path(salida_publica, "01_resumen_E05.csv"),
          row.names = FALSE, na = "")

campos_metricas <- c("total_conteos_genicos", "genes_con_conteos_positivos",
                     "genes_con_10_o_mas_conteos", "porcentaje_conteos_protein_coding")
lista_resumen <- lapply(campos_metricas, function(campo) {
  x <- metricas[[campo]]
  data.frame(metrica = campo, minimo = min(x, na.rm = TRUE),
             p25 = unname(quantile(x, 0.25, na.rm = TRUE)),
             mediana = median(x, na.rm = TRUE),
             media = mean(x, na.rm = TRUE),
             p75 = unname(quantile(x, 0.75, na.rm = TRUE)),
             maximo = max(x, na.rm = TRUE), stringsAsFactors = FALSE)
})
estadisticos <- do.call(rbind, lista_resumen)
write.csv(estadisticos, file.path(salida_publica, "02_estadisticos_QC.csv"),
          row.names = FALSE)

comparacion <- data.frame(
  indicador = c("Modelo de anotacion en bulk", "Tipo identificador GDC primario",
                "Tipo identificador GSE176078 en features",
                "Estado de lectura de referencia", "Detalle de error referencia",
                "Simbolos compartidos (ambos unicos)",
                "Ensembl IDs compartidos", "Uso del solapamiento"),
  resultado = c(modelo, "gene_id Ensembl con version",
                if (length(ref$ids_ensg)) "Ensembl IDs y simbolos" else
                  "simbolos/nombres (sin ENSG confirmado)",
                ref$estado, ref$error, as.character(length(comunes_simbolo)),
                if (is.na(id_compartidos)) "No calculable sin ENSG de referencia" else
                  as.character(id_compartidos),
                "Exploratorio; no fusionar matrices por simbolos sin resolver ambiguedades"),
  stringsAsFactors = FALSE)
write.csv(comparacion,
          file.path(salida_publica, "03_correspondencia_genes_GSE176078.csv"),
          row.names = FALSE)

metodo <- c(
  "# E/05 - Lectura y QC inicial de expresion bulk TCGA-BRCA", "",
  "**Estado:** resultados derivados de las 50 matrices STAR Counts descargadas y verificadas por E/04.",
  "Este script NO produce deconvolucion ni valida biologicamente proporciones celulares.", "",
  "## Entradas y procedencia", "",
  "- Cohorte: manifiesto FINAL E/03, corte GDC 2026-11-16, 50 archivos tumor primario.",
  "- Integridad: se exigen 50 MD5 y tamanos correctos; E/05 vuelve a calcular MD5.",
  "- Columnas GDC: `gene_id`, `gene_name`, `gene_type`, `unstranded`.",
  "- `unstranded` corresponde a conteos crudos. No se usan TPM ni FPKM.",
  "- Se excluyen de la matriz genica las cuatro filas tecnicas de STAR: `N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`.",
  "- No hay filtro de genes ni muestras en E/05; las metricas son descriptivas.", "",
  "## Resultados observados", "",
  sprintf("- Archivos bulk leidos: **%d**.", ncol(matriz)),
  sprintf("- Genes cuantificados por muestra: **%d**.", nrow(matriz)),
  sprintf("- Genes con conteos >0 en al menos una de las 50 muestras: **%d**.",
          sum(rowSums(matriz > 0) > 0)),
  sprintf("- Mediana de conteos genicos totales por muestra: **%.0f**.",
          median(metricas$total_conteos_genicos)),
  sprintf("- Mediana de genes detectados por muestra: **%.0f**.",
          median(metricas$genes_con_conteos_positivos)),
  sprintf("- Referencia GSE176078: **%s**.", ref$estado),
  sprintf("- Coincidencias exactas por simbolo, no ambiguas: **%d**.",
          length(comunes_simbolo)), "",
  "## Interpretacion y limitaciones", "",
  "La matriz local conserva los ID Ensembl con VERSION; no se sustituyen por simbolos.",
  "GSE176078 en 10x presenta nombres de genes; la coincidencia por simbolo es",
  "una evaluacion preliminar, no una armonizacion validada ni un mapeo uno-a-uno",
  "de IDs. Genes con nombres duplicados se excluyen de la comparacion por simbolo,",
  "pero NO de la matriz de conteos. Antes de B/C/E final se necesita resolver",
  "versiones, sinonimos, anotaciones y genes comunes de forma explicita.",
  "El total de conteos genicos NO representa el total de lecturas secuenciadas.",
  "No se elimina ninguna muestra por QC en esta etapa; revisar distribuciones",
  "y predefinir reglas antes de filtrar. No hay estimacion de celulas en bulk.", "",
  "## Ubicacion de resultados", "",
  "- Tablas agregadas para GitHub: `docs/E_Analisis_Bulk/05_expresion/`.",
  "- Graficos para GitHub: `figures/E_Analisis_Bulk/05_expresion/`.",
  "- Conteos, anotaciones y metricas individuales solo locales: `data/processed/TCGA_BRCA/05_expresion/`.",
  "- Entorno: `05_sessionInfo.txt`; codigo fuente: `scripts/E_Analisis_Bulk/05_leer_QC_expresion_TCGA_BRCA.R`.", "",
  "## Documentacion fuente", "",
  "- GDC mRNA Pipeline: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/"
)
writeLines(metodo, file.path(salida_publica, "04_metodologia_E05.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(salida_publica, "05_sessionInfo.txt"), useBytes = TRUE)

# Graficos descriptivos, sin IDs de paciente ni inferencia de biologica.
png(file.path(salida_fig, "01_conteos_por_muestra_hist.png"),
    width = 1200, height = 720, res = 140)
hist(log10(metricas$total_conteos_genicos + 1), breaks = "FD",
     col = "gray80", border = "gray40",
     main = "TCGA-BRCA (piloto): distribucion de conteos genicos",
     xlab = "log10(conteos unstranded por muestra + 1)",
     ylab = "Numero de muestras")
dev.off()

png(file.path(salida_fig, "02_genes_detectados_hist.png"),
    width = 1200, height = 720, res = 140)
hist(metricas$genes_con_conteos_positivos, breaks = "FD",
     col = "gray80", border = "gray40",
     main = "TCGA-BRCA (piloto): genes detectados por muestra",
     xlab = "Genes con conteo > 0", ylab = "Numero de muestras")
dev.off()

cat("\n=== E/05 | QC INICIAL DE EXPRESION BULK ===\n")
print(resumen, row.names = FALSE)
cat("\nDistribucion de metricas QC (sin identificadores de pacientes):\n")
print(estadisticos, row.names = FALSE)
cat("\nReferencia: ", ref$estado, "\n", sep = "")
if (nzchar(ref$error)) cat("Detalle referencia: ", ref$error, "\n", sep = "")
cat("\nArchivos publicables: ", salida_publica, " y ", salida_fig, "\n", sep = "")
cat("Matriz y metricas por muestra SOLO LOCALES: ", salida_local, "\n", sep = "")
cat("No se han estimado tipos celulares ni hecho deconvolucion.\n")
