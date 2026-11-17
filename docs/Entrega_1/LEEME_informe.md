# Entrega 1: informe reproducible en R Markdown

**Salida:** `Informe_Entrega_1.html` (HTML autónomo, sin Word).
**Fuente:** `Informe_Entrega_1.Rmd`, `estilo_informe.css` y `assets/flujo_metodologico.png`.
**Código de ejecución:** `scripts/00_render_informe_entrega1.R`.

1. Copiar las carpetas `docs/` y `scripts/` del paquete a la **raíz del repositorio** `proyecto-RNAseq-bulk`, conservando sus estructuras.
2. El repositorio debe contener las tablas existentes en `docs/A_Analisis_SingleCell/` y `docs/E_Analisis_Bulk/` de A/01–A/04 y E/01–E/05. Las tablas ya están publicadas en el repositorio; **no se repiten los análisis**.
3. En la consola de **R** de Positron (desde la raíz), comprobar dependencias:

   ```r
   requireNamespace("rmarkdown", quietly = TRUE)
   requireNamespace("knitr", quietly = TRUE)
   ```

   Si alguna respuesta es `FALSE`, revisar el entorno/`renv.lock` e instalar el paquete correspondiente. El lockfile del proyecto registra `rmarkdown` 2.32 y `knitr` 1.52, pero **no se ha probado una restauración limpia**.

4. Generar HTML:

   ```r
   source("scripts/00_render_informe_entrega1.R")
   ```

5. Abrir `docs/Entrega_1/Informe_Entrega_1.html`. Si se requiere entregar como PDF, usar la opción *Imprimir → Guardar como PDF* del navegador, verificando **tres páginas A4**. Esto es opcional; el informe solicitado aquí no genera Word.
6. Publicar las fuentes y, tras comprobarlo, el HTML en GitHub:

   ```powershell
   git add docs/Entrega_1/ scripts/00_render_informe_entrega1.R
   git diff --cached --check
   git status
   git commit -m "Entrega 1: informe reproducible R Markdown y flujo metodologico"
   git push origin main
   ```

## Trazabilidad

- `docs/A_Analisis_SingleCell/01_objeto_seurat/resumen.csv`: 26 donantes, 100.064 células, 29.733 genes.
- `docs/A_Analisis_SingleCell/04_filtrado_qc/resumen_conservacion.csv`: decisión de conservación/QC inicial.
- `docs/E_Analisis_Bulk/01_cohorte/01_resumen_TCGA_BRCA.csv`: 1.095 pacientes, 1.111 archivos, elegibilidad no FFPE.
- `docs/E_Analisis_Bulk/02_clinica/01_resumen_clinico.csv`: disponibilidad clínica.
- `docs/E_Analisis_Bulk/03_piloto/`: 50 pacientes, estadio y manifiesto con tamaño/MD5.
- `docs/E_Analisis_Bulk/04_descarga/`: integridad 50/50 archivos.
- `docs/E_Analisis_Bulk/05_expresion/`: 50 matrices leídas, 60.660 genes, 20.288 símbolos compartidos, estadísticas descriptivas.

El informe **se detiene ante inconsistencias** del corte fijo del 16-11-2026 (falta de archivos, discordancia de paciente/manifest/MD5, conteos críticos). No descarga matrices, no estima composición celular y no valida el rendimiento de MuSiC/NNLS.

### Limitaciones de reproducibilidad

- Los datos originales y matrices en `data/raw/` o `data/processed/` permanecen fuera del repositorio. El informe se regenera desde **CSV agregados ya versionados**.
- `renv.lock` está versionado, pero su restauración desde un clon limpio y la portabilidad de BPCells aún deben probarse.
- El HTML puede tener pequeñas diferencias de paginación entre navegadores y fuentes; revisar la exportación visual antes de enviar.
- El diagrama es una **representación gráfica del procedimiento previsto**: la etapa B–D todavía no está ejecutada.
