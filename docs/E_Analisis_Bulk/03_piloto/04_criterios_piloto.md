# Piloto TCGA-BRCA: cierre de la seleccion (E03)

Fecha de los metadatos originales: 2026-11-16.
E01: 50 pacientes seleccionados aleatoriamente mediante set.seed(20261116),
despues de una muestra primaria por Case ID y exclusión de FFPE.
E03 NO selecciona por estadio: solamente comprueba su disponibilidad.
Se conservan exactamente los mismos 50 File ID de E01 (sin reemplazar casos).
Fuente clinica: E02, solo diagnostico de enfermedad primaria.
Tamanio estimado: 212.050 MB decimales (suma size del manifest).

## Limitaciones
El piloto es de factibilidad tecnica y NO tiene potencia para comparar estadios.
No confundir 50 muestras con validacion biológica ni con matrices descargadas.
La variable Unknown en preservacion no equivale a OCT confirmado.
El estadio IV queda pobremente representado; no realizar inferencias clinicas fuertes.
El GDC puede actualizar archivos: conservar manifest, fecha y MD5.

## Paso siguiente
Instalar gdc-client (Windows), descargar -m 03_manifest_piloto_50_FINAL.txt
en data/raw/TCGA_BRCA/expresion_piloto (excluido de Git) y ejecutar E04.
Verificar genes compartidos y escala de cuantificacion antes de deconvolucion.
Se registra MD5 de entradas GDC en 05_md5_insumos.csv.
