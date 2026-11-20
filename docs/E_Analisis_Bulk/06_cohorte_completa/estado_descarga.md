# TCGA-BRCA — Estado de cohorte completa

Fecha: 2026-11-19.

## Cohorte seleccionada
- 1093 pacientes no FFPE.
- 679 muestras OCT y 414 de preservación Unknown.
- Una muestra tumoral primaria por paciente.
- Manifiesto de 1093 archivos; 4.638 GB previstos.
- Los 50 pacientes del piloto están incluidos.

## Estado de datos
- Cohorte y manifiesto: preparados y validados.
- Piloto original de 50 archivos: descarga y QC realizados previamente.
- Cohorte ampliada: 1093 archivos descargados en Windows.
- Integridad: 1093/1093 archivos con MD5 correcto.
- Evidencia: `integridad_1093.csv` y script PowerShell versionado.
- Transferencia a Kutral y QC de las 1093 matrices: pendientes.
- No existe todavía un QC de las 1093 matrices completas.
- No se ha realizado deconvolución.

## Incidencia técnica
El 19-11-2026 se comprobaron fallos de conexión al
descargar archivos desde GDC, tanto en un nodo de
cómputo de Kutral como mediante gdc-client en Windows.

El portal GDC mostró un aviso de intermitencia de
herramientas de análisis, sin confirmación de que
esa incidencia sea la causa de los errores de descarga.

Se conserva el manifiesto y no se modifican los criterios
de selección ni los resultados del piloto.
