# E/04 — Protocolo reproducible de descarga e integridad (TCGA-BRCA)

## Contexto y corte de datos

- Fuente: GDC, archivos públicos `TCGA-BRCA` / RNA-seq / `STAR - Counts` / tumor primario.
- Exportación de metadatos GDC: **16-11-2026**.
- Manifiesto *final* generado por `E/03`: `docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt`.
- Cohorte piloto: 50 archivos, una muestra por paciente; 212.050.127 bytes esperados.
- Herramienta efectivamente usada por la investigadora: **GDC Data Transfer Tool 2.3**, ejecutada en PowerShell sobre Windows.
- La descarga original fue realizada manualmente con `gdc-client.exe` siguiendo las instrucciones en conversación. El script `04_descargar_bulk_TCGA_BRCA.ps1` se añade **después** como procedimiento automatizado para nuevas ejecuciones; no afirmar que ese script en sí ya fue ejecutado.

## Procedimiento

1. Descargar y descomprimir `gdc-client` desde https://gdc.cancer.gov/access-data/gdc-data-transfer-tool. Comprobar `--version` (en el equipo original: `2.3`).
2. Asegurar que el manifiesto E/03 existe y tiene 50 filas de datos, IDs únicos y campos `id`, `filename`, `md5`, `size`.
3. Desde PowerShell en la raíz del repositorio, ejecutar:

   ```powershell
   & "scripts/E_Analisis_Bulk/04_descargar_bulk_TCGA_BRCA.ps1"
   ```

   Si no localiza automáticamente el ejecutable, indicar `-GdcExe "C:\ruta\gdc-client.exe"`.
   El script crea el directorio `data/raw/TCGA_BRCA/expresion_piloto/` antes de descargar, guarda SHA256 del manifiesto y ejecutable (más versión y fecha UTC) en un registro JSON bajo `docs/E_Analisis_Bulk/04_descarga/`.
4. Ejecutar **en R**:

   ```r
   source('scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R')
   ```

5. Consultar `01_integridad_archivos.csv` y `02_resumen_descarga.csv`. Declarar descarga verificada **solo si** hay 50 presentes, 50 MD5/tamaños válidos y cero errores.

## Resultado observado en el equipo original

El usuario ejecutó `04_verificar_descarga_TCGA_BRCA.R` después de descargar y obtuvo:

| Indicador | Valor observado |
|---|---:|
| Archivos esperados | 50 |
| Archivos presentes | 50 |
| Archivos íntegros por MD5 y tamaño | 50 |
| Ausentes o con error | 0 |
| Bytes esperados | 212050127 |
| Estado | Verificada |

El encabezado de un archivo confirmó `# gene-model: GENCODE v36` y columnas `gene_id`, `gene_name`, `gene_type`, `unstranded`, `stranded_first`, `stranded_second`, `tpm_unstranded`, `fpkm_unstranded`, `fpkm_uq_unstranded`. Esto permite planificar E/05, **no equivale a QC cuantitativo de expresión**.

## Separación de responsabilidades

- PowerShell `.ps1`: **obtención** de los datos crudos.
- R E/04: **verificación y reporte** de integridad.
- R E/05: **lectura, filtrado y QC** de las matrices; todavía pendiente.
- Los datos crudos `data/raw/` quedan fuera de GitHub por `.gitignore`.
- No se publican rutas locales con nombre de usuario ni tokens de acceso.

## Limitaciones reproducibles

- El servidor GDC y su contenido pueden cambiar; la exportación y su fecha fijan el corte original. Para reproducir exactamente los mismos archivos, usar el manifiesto versionado, no una búsqueda nueva.
- El registro JSON generado por el nuevo script documenta la ejecución futura, pero no sustituye los MD5 verificados por R.
- Para reproducibilidad de extremo a extremo se necesita además el acceso a la referencia GSE176078, ejecutar A y los pasos B–D, y probar `renv::restore()` en entorno limpio; esto último no se ha verificado todavía.
