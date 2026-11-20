# Entrega 1 — informe reproducible en R Markdown

- **Fuente:** `Informe_Entrega_1.Rmd` y CSV versionados en `docs/A_Analisis_SingleCell/` y `docs/E_Analisis_Bulk/`.
- **Script generador:** `scripts/00_render_informe_entrega1.R`.
- **Salida:** `Informe_Entrega_1.html`.
- **Figura:** `assets/flujo_metodologico.png`, incrustada en Base64.
- **Estilos:** `estilo_informe.css`, integrados en el HTML.

El HTML generado es autónomo. Se conserva el PNG original
para poder regenerar el documento.

## Regeneración en Kutral

Desde la raíz del proyecto, ejecutar mediante Slurm:

    sbatch scripts/00_render_informe_entrega1.sbatch

Se utiliza R 4.5.1, rmarkdown, knitr, base64enc y Pandoc 3.14.
El informe se genera desde resultados y manifiestos versionados.
No descarga datos ni realiza deconvolución.

## Revisión visual

Abrir el HTML en el navegador y comprobar la impresión A4.
Revisar que el documento ocupe como máximo tres páginas
y que las tablas e imágenes no estén cortadas.

## Archivos a publicar

```powershell
git add docs/Entrega_1/Informe_Entrega_1.Rmd docs/Entrega_1/Informe_Entrega_1.html
git add docs/Entrega_1/LEEME_informe.md docs/Entrega_1/assets/flujo_metodologico.png
git add docs/Entrega_1/estilo_informe.css scripts/00_render_informe_entrega1.R
git diff --cached --check
git status
```

Los datos originales (`data/raw/`) y objetos grandes (`data/processed/`) permanecen fuera de GitHub. Los scripts y las tablas publicadas proporcionan trazabilidad de los resultados ya obtenidos. La restauración limpia de `renv` y la portabilidad de BPCells no se han validado todavía.
