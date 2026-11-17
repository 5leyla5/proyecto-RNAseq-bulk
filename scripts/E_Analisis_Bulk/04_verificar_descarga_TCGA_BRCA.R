# =============================================================================
# Bioinformatica II | TCGA-BRCA / E04: integridad de descarga STAR Counts
# Se ejecuta DESPUES de usar gdc-client. No descarga, altera ni carga 50 matrices.
# source('scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R')
# =============================================================================
stopifnot(file.exists('README.md'))
manifest <- file.path('docs','E_Analisis_Bulk','03_piloto','03_manifest_piloto_50_FINAL.txt')
base <- file.path('data','raw','TCGA_BRCA','expresion_piloto')
salida <- file.path('docs','E_Analisis_Bulk','04_descarga')
if (!file.exists(manifest)) stop('Primero ejecutar E03: falta el manifest FINAL.')
m <- read.delim(manifest,sep='\t',stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(m)==50L, !anyDuplicated(m$id))
dir.create(salida,recursive=TRUE,showWarnings=FALSE)
# gdc-client guarda cada archivo bajo su UUID: -d base/<uuid>/<filename>
rutas <- file.path(base,m$id,m$filename)
existe <- file.exists(rutas)
bytes_locales <- rep(NA_real_,length(rutas))
md5_locales <- rep(NA_character_,length(rutas))
if (any(existe)) {
  bytes_locales[existe] <- file.info(rutas[existe])$size
  md5_locales[existe] <- unname(tools::md5sum(rutas[existe]))
}
valido <- existe & !is.na(bytes_locales) &
  bytes_locales == as.numeric(m$size) &
  !is.na(md5_locales) & tolower(md5_locales)==tolower(m$md5)
resultado <- data.frame(file_id=m$id,nombre=m$filename,
  existe=existe,size_esperado=m$size,size_local=bytes_locales,
  md5_esperado=m$md5,md5_local=md5_locales,integro=valido,
  stringsAsFactors=FALSE)
write.csv(resultado,file.path(salida,'01_integridad_archivos.csv'),row.names=FALSE,na='')
resumen <- data.frame(indicador=c('Archivos esperados','Archivos presentes',
  'Archivos integro MD5 y tamano','Ausentes o con error','Bytes esperados',
  'Estado de descarga completa'),
  valor=c(nrow(m),sum(existe),sum(valido),sum(!valido),sum(as.numeric(m$size)),
          if(all(valido)) 'Verificada' else 'Pendiente / revisar'),stringsAsFactors=FALSE)
write.csv(resumen,file.path(salida,'02_resumen_descarga.csv'),row.names=FALSE)
cat('\n=== E04 | INTEGRIDAD DE DESCARGA ===\n'); print(resumen,row.names=FALSE)
if (all(valido)) {
  # Inspeccionar encabezado de solo UNA matriz para orientar el preprocesamiento.
  muestra <- readLines(rutas[[1]],n=7L,warn=FALSE)
  cat('\nPrimeras lineas del primer archivo STAR Counts:\n')
  cat(paste(muestra,collapse='\n'),'\n')
  cat('\nArchivos completos. Siguiente etapa: lectura con gene_id / unstranded,\n',
      'correspondencia de identificadores geneticos y QC de expresion.\n',sep='')
} else {
  cat('\nNo declarar descarga terminada. Reanuda gdc-client con el mismo manifest.\n')
}
