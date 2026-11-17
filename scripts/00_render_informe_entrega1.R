# Genera un HTML reproducible de la Entrega 1 a partir de CSV versionados.
# Ejecutar en la consola R de Positron, desde la raiz del repositorio:
# source("scripts/00_render_informe_entrega1.R")
# Requiere rmarkdown y knitr; versiones registradas en renv.lock.

if (!file.exists("README.md") || !dir.exists("docs/E_Analisis_Bulk")) {
  stop("Abra la raiz del repositorio Proyecto RNA-seq bulk antes de generar el informe.")
}
paquetes <- c("rmarkdown", "knitr")
faltantes <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes)) {
  stop("Faltan paquetes: ", paste(faltantes, collapse = ", "),
       ". Instálelos en el entorno R del proyecto y registre sus versiones en renv.")
}
entrada <- "docs/Entrega_1/Informe_Entrega_1.Rmd"
if (!file.exists(entrada)) stop("Falta archivo: ", entrada)
if (!file.exists("docs/Entrega_1/assets/flujo_metodologico.png")) {
  stop("Falta la imagen del flujo metodologico.")
}
raiz <- normalizePath(getwd(), winslash = "/")
resultado <- rmarkdown::render(
  input = entrada,
  output_format = "html_document",
  output_file = "Informe_Entrega_1.html",
  output_dir = file.path(raiz, "docs", "Entrega_1"),
  knit_root_dir = raiz,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
cat("\nInforme HTML generado desde los CSV del proyecto:\n", normalizePath(resultado, winslash = "/"), "\n", sep = "")
cat("Revisar visualmente el HTML y, si se imprime, comprobar paginacion A4. No se genero Word.\n")
