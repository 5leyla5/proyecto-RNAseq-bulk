# Entrega 1 - renderizado reproducible HTML (sin procesar PNG como texto).
# Ejecutar desde la raiz del proyecto en la consola R de Positron:
# source("scripts/00_render_informe_entrega1.R", encoding = "UTF-8")
# El HTML usa assets/flujo_metodologico.png como archivo vinculado.
# Mantener el HTML y su carpeta assets/ juntos en GitHub / al compartir.

if (!file.exists("README.md") || !dir.exists("docs/E_Analisis_Bulk")) {
  stop("Abre la raiz del repositorio Proyecto RNA-seq bulk en Positron.")
}

paquetes <- c("rmarkdown", "knitr")
faltantes <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes)) {
  stop("Paquetes R faltantes: ", paste(faltantes, collapse = ", "))
}

raiz <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
entrada <- file.path(raiz, "docs", "Entrega_1", "Informe_Entrega_1.Rmd")
imagen <- file.path(raiz, "docs", "Entrega_1", "assets", "flujo_metodologico.png")
css <- file.path(raiz, "docs", "Entrega_1", "estilo_informe.css")
salida <- file.path(raiz, "docs", "Entrega_1", "Informe_Entrega_1.html")

for (p in c(entrada, imagen, css)) {
  if (!file.exists(p)) stop("Falta archivo de entrada: ", p)
}
# Validacion de cabecera PNG; no se interpreta ni lee el archivo como texto.
con <- file(imagen, open = "rb")
firma <- readBin(con, what = "raw", n = 8L)
close(con)
firma_png <- as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))
if (!identical(firma, firma_png)) stop("La imagen no tiene cabecera PNG valida.")

tmp <- tempfile("entrega1_fragmento_")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
cat("[1/2] Generando texto y tablas HTML desde los CSV...\n")
fragmento <- rmarkdown::render(
  input = entrada,
  output_format = rmarkdown::html_fragment(),
  output_file = "fragmento.html",
  output_dir = tmp,
  knit_root_dir = raiz,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

texto <- paste(readLines(fragmento, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
marca <- "<!-- INSERTAR_FIGURA_FLUJO_POST_PANDOC -->"
if (!grepl(marca, texto, fixed = TRUE)) {
  stop("No esta el marcador de imagen en el Rmd. Debe contener INSERTAR_FIGURA_FLUJO_POST_PANDOC.")
}

# No crear dataURI ni incrustar la imagen como una cadena Base64 gigante.
# Esta ruta se resuelve desde el HTML en docs/Entrega_1/.
cat("[2/2] Enlazando assets/flujo_metodologico.png e insertando estilos...\n")
etiqueta <- paste0(
  '<img src="assets/flujo_metodologico.png" ',
  'alt="Flujo metodologico: GSE176078, referencia, pseudobulk, NNLS y MuSiC, validacion y TCGA-BRCA" ',
  'style="max-height:9.3cm;max-width:100%;height:auto;" />'
)
texto <- sub(marca, etiqueta, texto, fixed = TRUE)
estilos <- paste(readLines(css, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (nchar(estilos) < 400L || !grepl("font-family", estilos, fixed = TRUE)) {
  stop("La hoja de estilos esta vacia o no es valida: ", css)
}
if (grepl("<!DOCTYPE|<html|<head|<body", texto, ignore.case = TRUE)) {
  stop("El fragmento Pandoc incluye etiquetas de documento completo: revisa la salida html_fragment.")
}
# Proteccion visual: estilos base integrados, aun si un navegador no aplica
# el CSS externo. Se fuerzan tipografias y contenedores en pantalla.
estilo_respaldo <- paste(c(
  'html, body {font-family: "Segoe UI", Arial, sans-serif !important; color:#1b2635 !important; background:#f3f6fb !important; font-size:10.5pt !important; line-height:1.36 !important;}',
  'body {margin:0 !important;}',
  '.main-container {max-width:820px !important; margin:24px auto !important; padding:0 !important;}',
  '#header {background:white !important; padding:18px 24px 8px !important; border-top:5px solid #384c81 !important; margin-bottom:12px !important; border-radius:8px !important;}',
  '#header h1 {font-family:"Segoe UI",Arial,sans-serif !important; font-size:21pt !important; font-weight:700 !important; line-height:1.17 !important; margin:0 0 6px !important;}',
  '#header .author,#header .date {font-size:9pt !important; color:#586273 !important;}',
  '.report-page {font-family:"Segoe UI",Arial,sans-serif !important; background:#fff !important; padding:18px 24px 13px !important; margin:0 0 15px !important; border:1px solid #dbe0e8 !important; border-radius:8px !important;}',
  '.report-page h2 {font-family:"Segoe UI",Arial,sans-serif !important; font-size:15pt !important; margin:5px 0 10px !important; padding-bottom:5px !important; border-bottom:2px solid #c7d0df !important;}',
  '.report-page h3 {font-family:"Segoe UI",Arial,sans-serif !important; font-size:11.5pt !important; color:#31476c !important; margin:12px 0 6px !important;}',
  '.report-page p {margin:6px 0 9px !important;}',
  '.report-page table {border-collapse:collapse !important; font-size:8.7pt !important; width:100% !important;}',
  '.report-page td,.report-page th {border:1px solid #cad3de !important; padding:4px 6px !important;}',
  '.report-page th {background:#e9edf5 !important;}',
  '.figure {text-align:center !important;}',
  '.figure img {max-height:9.3cm !important; max-width:100% !important; width:auto !important;}',
  '@media print {',
  ' @page {size:A4 portrait; margin:11mm 13mm 10mm;}',
  ' html,body {background:white !important; font-size:9.5pt !important; line-height:1.25 !important;}',
  ' .main-container {max-width:none !important; margin:0 !important;}',
  ' #header {border:0 !important; border-top:3px solid #344b7c !important; padding:0 0 5px !important; margin:0 0 5px !important; border-radius:0 !important;}',
  ' #header h1 {font-size:16pt !important;}',
  ' .report-page {background:white !important; padding:0 !important; margin:0 !important; border:none !important; border-radius:0 !important; break-after:page; page-break-after:always;}',
  ' .report-page:last-of-type {break-after:auto; page-break-after:auto;}',
  ' .report-page h2 {font-size:13pt !important; margin:4px 0 6px !important;}',
  ' .report-page h3 {font-size:10.5pt !important; margin:8px 0 4px !important;}',
  ' .report-page p {margin:4px 0 7px !important;}',
  ' .report-page td,.report-page th {padding:3px 5px !important;}',
  ' .figure img {max-height:9.0cm !important;}',
  ' .report-page h2,.report-page h3,.figure,table {break-inside:avoid;page-break-inside:avoid;}',
  '}'
), collapse = "\n")

titulo <- paste0(
  "Diseño y factibilidad de un flujo reproducible para ",
  "deconvolución celular de RNA-seq bulk en cáncer de mama"
)
cabecera <- paste0(
  '<header id="header">',
  '<h1>', titulo, '</h1>',
  '<p class="author">Leyla Victoria Saavedra Quinteros · ',
  'Magíster en Informática Médica · Universidad de Chile</p>',
  '<p class="date">Bioinformática II · Entrega parcial 1 · 19 de noviembre de 2026</p>',
  '</header>'
)
html <- paste0(
  '<!DOCTYPE html>\n<html lang="es">\n<head>\n',
  '<meta charset="UTF-8" />\n',
  '<meta name="viewport" content="width=device-width, initial-scale=1" />\n',
  '<title>', titulo, '</title>\n',
  '<link rel="stylesheet" href="estilo_informe.css" />\n',
  '<style type="text/css">\n', estilos, '\n', estilo_respaldo, '\n</style>\n',
  '</head>\n<body>\n<div class="main-container">\n',
  cabecera, '\n', texto, '\n',
  '</div>\n</body>\n</html>\n'
)
if (!grepl('<style type="text/css">', html, fixed = TRUE) ||
    !grepl('font-family', html, fixed = TRUE) ||
    !grepl('href="estilo_informe.css"', html, fixed = TRUE)) {
  stop("Fallo la inclusion de los estilos del informe.")
}
writeLines(enc2utf8(html), con = salida, useBytes = TRUE)
if (!file.exists(salida) || file.info(salida)$size < 1000) {
  stop("No se genero correctamente el HTML.")
}
# No leer de nuevo el archivo completo mediante grep()/readLines():
# esa comprobacion originaba el error de codificacion en la version anterior.
secciones <- gregexpr('class="[^"]*report-page', texto, perl = TRUE)[[1L]]
n_paginas <- if (identical(secciones[1L], -1L)) 0L else length(secciones)
if (n_paginas != 3L) {
  warning("Se encontraron ", n_paginas, " secciones de informe; revisar paginacion A4.")
}
unlink(tmp, recursive = TRUE, force = TRUE)
cat("\nInforme HTML generado: ", salida,
    "\nFigura enlazada correctamente: docs/Entrega_1/assets/flujo_metodologico.png",
    "\nMantener HTML, estilo_informe.css y assets/ juntos. No se genero Word.",
    "\nComprobar visualmente 3 paginas A4 antes del commit.\n", sep = "")
