# Entrega 1 — informe reproducible en R Markdown

- **Fuente:** `Informe_Entrega_1.Rmd` y CSV versionados en `docs/A_Analisis_SingleCell/` y `docs/E_Analisis_Bulk/`.
- **Script generador:** `scripts/00_render_informe_entrega1.R`.
- **Salida:** `Informe_Entrega_1.html`.
- **Figura:** `assets/flujo_metodologico.png` (**archivo vinculado**, no incrustado en Base64).
- **Estilos:** `estilo_informe.css`, integrados en el HTML por el script.

La imagen se carga mediante una ruta relativa para evitar errores de codificación al procesar un PNG. Por ello, **el HTML no es un archivo único completamente autónomo**: debe conservarse junto a `assets/flujo_metodologico.png` dentro de `docs/Entrega_1/` al compartir o publicar.

## Regeneración

Abrir Positron en la raíz del repositorio y ejecutar en la consola R:

```r
source("scripts/00_render_informe_entrega1.R", encoding = "UTF-8")
```

Se requieren los paquetes R `rmarkdown` y `knitr` y los CSV versionados ya existentes. El informe no descarga datos, no recalcula E/01–E/05, no aplica MuSiC/NNLS ni genera Word. El script genera primero un fragmento HTML desde las tablas, después inserta el enlace relativo a la figura y los estilos CSS.

Abrir desde PowerShell:

```powershell
Start-Process "docs/Entrega_1/Informe_Entrega_1.html"
```

Comprobar que se vea la figura en la sección 2 y revisar la vista de impresión **A4, máximo tres páginas**, desactivando encabezados y pies del navegador si corresponde.

## Archivos a publicar

```powershell
git add docs/Entrega_1/Informe_Entrega_1.Rmd docs/Entrega_1/Informe_Entrega_1.html
git add docs/Entrega_1/LEEME_informe.md docs/Entrega_1/assets/flujo_metodologico.png
git add docs/Entrega_1/estilo_informe.css scripts/00_render_informe_entrega1.R
git diff --cached --check
git status
```

Los datos originales (`data/raw/`) y objetos grandes (`data/processed/`) permanecen fuera de GitHub. Los scripts y las tablas publicadas proporcionan trazabilidad de los resultados ya obtenidos. La restauración limpia de `renv` y la portabilidad de BPCells no se han validado todavía.
