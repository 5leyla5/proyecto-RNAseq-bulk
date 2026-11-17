# =============================================================================
# Bioinformatica II | Proyecto RNA-seq bulk | TCGA-BRCA / E03
# Cierra piloto de 50 pacientes a partir de E01 y E02; NO descarga matrices.
# Ejecutar desde raiz del proyecto en consola R de Positron:
# source('scripts/E_Analisis_Bulk/03_cerrar_piloto_TCGA_BRCA.R')
# =============================================================================
stopifnot(file.exists('README.md'), dir.exists('scripts/E_Analisis_Bulk'))
raiz <- file.path('docs','E_Analisis_Bulk')
r_piloto <- file.path(raiz,'01_cohorte','04_piloto_50_provisional.csv')
r_manifest <- file.path(raiz,'01_cohorte','05_manifest_piloto_50.txt')
r_clinico <- file.path('data','processed','TCGA_BRCA','clinica_por_paciente.csv')
for (archivo in c(r_piloto, r_manifest, r_clinico)) {
  if (!file.exists(archivo)) stop('Falta ', archivo, '. Ejecuta primero E01 y E02.')
}
leer_csv <- function(r) read.csv(r, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('NA',''))
p <- leer_csv(r_piloto)
m <- read.delim(r_manifest, sep='\t', check.names=FALSE, stringsAsFactors=FALSE)
clin <- leer_csv(r_clinico)
stopifnot(nrow(p)==50L, !anyDuplicated(p[['Case ID']]), !anyDuplicated(p[['File ID']]),
          !anyDuplicated(clin$case_id), nrow(m)==50L, !anyDuplicated(m$id))
if (!all(p[['Project ID']]=='TCGA-BRCA' & p[['Tissue Type']]=='Tumor' &
         p[['Tumor Descriptor']]=='Primary' & p[['Preservation Method']]!='FFPE')) {
  stop('El piloto no cumple los criterios establecidos en E01.')
}
if (!setequal(p[['File ID']],m$id)) stop('Piloto/manifest E01 no coinciden en File ID.')
ixm <- match(p[['File ID']],m$id)
if (!identical(as.character(p[['File Name']]),as.character(m$filename[ixm])) ||
    !identical(as.character(p$md5),as.character(m$md5[ixm])) ||
    !all(as.numeric(p$size_bytes)==as.numeric(m$size[ixm]))) {
  stop('Discordancia entre catalogo del piloto y manifest original.')
}
ixc <- match(p[['Case ID']],clin$paciente_tcga)
if (anyNA(ixc)) stop('Faltan pacientes del piloto en la tabla clinica E02.')
cp <- clin[ixc,,drop=FALSE]
if (!all(cp$estado_estadio=='Disponible') || anyNA(cp$estadio_agrupado) ||
    any(!cp$estadio_agrupado %in% c('I','II','III','IV'))) {
  stop('Existen pacientes sin estadio AJCC utilizable. Mantener E01 provisional y revisar sin reemplazos ad hoc.')
}
if (!all(p[['Case ID']]==cp$paciente_tcga)) stop('Cruce clinico incorrecto.')
# No se usa el estadio para seleccionar, reemplazar ni equilibrar pacientes.
# Se conserva exactamente la muestra aleatoria seleccionada ANTES de ver la clinica.
datos <- data.frame(case_id=cp$case_id, paciente_tcga=cp$paciente_tcga,
                    sample_id=p[['Sample ID']], file_id=p[['File ID']],
                    estadio_ajcc=cp$estadio_ajcc,
                    estadio_agrupado=cp$estadio_agrupado,
                    edad_anos=cp$edad_anos, preservacion=p[['Preservation Method']],
                    stringsAsFactors=FALSE)
out <- file.path(raiz,'03_piloto')
local <- file.path('data','processed','TCGA_BRCA')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(local,recursive=TRUE,showWarnings=FALSE)
write.csv(datos,file.path(local,'piloto_50_clinica.csv'),row.names=FALSE,na='')
conteos <- table(factor(datos$estadio_agrupado,levels=c('I','II','III','IV')))
write.csv(data.frame(estadio=names(conteos),pacientes=as.integer(conteos)),
          file.path(out,'02_estadios_piloto.csv'),row.names=FALSE)
mb <- sum(as.numeric(m$size))/1e6
resumen <- data.frame(indicador=c('Pacientes unicos del piloto',
                                  'Muestras tumor primario',
                                  'Pacientes con estadio AJCC utilizable',
                                  'Archivos OCT','Archivos Unknown','Archivos FFPE',
                                  'Tamano estimado de descarga (MB decimales)',
                                  'Archivos de expresion descargados por E03'),
                      valor=c(50L,50L,sum(cp$estado_estadio=='Disponible'),
                              sum(p[['Preservation Method']]=='OCT'),
                              sum(p[['Preservation Method']]=='Unknown'),
                              sum(p[['Preservation Method']]=='FFPE'),
                              round(mb,3),NA),stringsAsFactors=FALSE)
write.csv(resumen,file.path(out,'01_resumen_piloto.csv'),row.names=FALSE,na='')
# Manifest listo para cliente GDC, en el ORDEN del piloto E01.
manifest_final <- m[ixm,c('id','filename','md5','size','state')]
write.table(manifest_final,file.path(out,'03_manifest_piloto_50_FINAL.txt'),
            sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
notas <- c('# Piloto TCGA-BRCA: cierre de la seleccion (E03)','',
 'Fecha de los metadatos originales: 2026-11-16.',
 'E01: 50 pacientes seleccionados aleatoriamente mediante set.seed(20261116),',
 'despues de una muestra primaria por Case ID y exclusión de FFPE.',
 'E03 NO selecciona por estadio: solamente comprueba su disponibilidad.',
 'Se conservan exactamente los mismos 50 File ID de E01 (sin reemplazar casos).',
 'Fuente clinica: E02, solo diagnostico de enfermedad primaria.',
 sprintf('Tamanio estimado: %.3f MB decimales (suma size del manifest).',mb),'',
 '## Limitaciones',
 'El piloto es de factibilidad tecnica y NO tiene potencia para comparar estadios.',
 'No confundir 50 muestras con validacion biológica ni con matrices descargadas.',
 'La variable Unknown en preservacion no equivale a OCT confirmado.',
 'El estadio IV queda pobremente representado; no realizar inferencias clinicas fuertes.',
 'El GDC puede actualizar archivos: conservar manifest, fecha y MD5.','',
 '## Paso siguiente',
 'Instalar gdc-client (Windows), descargar -m 03_manifest_piloto_50_FINAL.txt',
 'en data/raw/TCGA_BRCA/expresion_piloto (excluido de Git) y ejecutar E04.',
 'Verificar genes compartidos y escala de cuantificacion antes de deconvolucion.',
 'Se registra MD5 de entradas GDC en 05_md5_insumos.csv.')
writeLines(notas,file.path(out,'04_criterios_piloto.md'),useBytes=TRUE)
# Huellas digitales de las ENTRADAS locales (no se versionan datos crudos).
# Conservar estos MD5 junto con el manifiesto permite auditar la instantanea GDC.
entradas <- c(file.path('data','raw','TCGA_BRCA','gdc_sample_sheet.tsv'),
              file.path('data','raw','TCGA_BRCA','gdc_manifest.txt'),
              file.path('data','raw','TCGA_BRCA','cohort_TCGA_BRCA_LSQ.2026-11-16.tsv'),
              file.path('data','raw','TCGA_BRCA','clinical','clinical.tsv'))
if (!all(file.exists(entradas))) stop('Faltan entradas para registrar MD5: ',
  paste(entradas[!file.exists(entradas)],collapse=', '))
registro_md5 <- data.frame(archivo=entradas,
                           md5=unname(tools::md5sum(entradas)),
                           stringsAsFactors=FALSE)
write.csv(registro_md5,file.path(out,'05_md5_insumos.csv'),row.names=FALSE)
cat('\n=== E03 | PILOTO CERRADO, NO DESCARGADO ===\n')
print(resumen,row.names=FALSE)
cat('\nEstadios del piloto:\n'); print(conteos)
cat('\nManifest: ',file.path(out,'03_manifest_piloto_50_FINAL.txt'),'\n',sep='')
cat('Detalle por paciente (SOLO LOCAL): ',file.path(local,'piloto_50_clinica.csv'),'\n',sep='')
