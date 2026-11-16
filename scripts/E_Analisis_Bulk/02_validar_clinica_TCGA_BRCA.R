# =============================================================================
# Bioinformatica II | Proyecto RNA-seq bulk
# Objetivo E / 02 - Validacion clinica y estadios AJCC de TCGA-BRCA
# Entradas: data/raw/TCGA_BRCA/clinical/clinical.tsv
#           data/raw/TCGA_BRCA/cohort_TCGA_BRCA_LSQ.2026-11-16.tsv
# Salidas publicables: docs/E_Analisis_Bulk/02_clinica/*.csv y *.md
# Salida local excluida de Git: data/processed/TCGA_BRCA/clinica_por_paciente.csv
# Ejecutar DESDE LA RAIZ del proyecto en la consola R de Positron:
#   source('scripts/E_Analisis_Bulk/02_validar_clinica_TCGA_BRCA.R')
# Solo usa R base. No descarga datos ni modifica los datos de GSE176078.
# =============================================================================

if (!file.exists('README.md') || !dir.exists('scripts/E_Analisis_Bulk')) {
  stop('Ejecuta desde la raiz del proyecto Proyecto RNA-seq bulk.')
}

raiz_datos <- file.path('data', 'raw', 'TCGA_BRCA')
ruta_clinica <- file.path(raiz_datos, 'clinical', 'clinical.tsv')
ruta_cohorte <- file.path(raiz_datos, 'cohort_TCGA_BRCA_LSQ.2026-11-16.tsv')
salidas <- file.path('docs', 'E_Analisis_Bulk', '02_clinica')
salida_local <- file.path('data', 'processed', 'TCGA_BRCA')

for (p in c(ruta_clinica, ruta_cohorte)) {
  if (!file.exists(p)) stop('No encuentro: ', p,
                             '\nComprueba la carpeta y la descompresion de GDC.')
}
dir.create(salidas, recursive = TRUE, showWarnings = FALSE)
dir.create(salida_local, recursive = TRUE, showWarnings = FALSE)

leer_tsv <- function(p) {
  utils::read.delim(p, sep = '\t', check.names = FALSE, quote = '',
                    stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
}
cl <- leer_tsv(ruta_clinica)
coh <- leer_tsv(ruta_cohorte)

campos <- c('project.project_id', 'cases.case_id', 'cases.submitter_id',
            'diagnoses.diagnosis_is_primary_disease',
            'diagnoses.ajcc_pathologic_stage', 'demographic.age_at_index')
if (!all(campos %in% names(cl))) {
  stop('Faltan columnas: ', paste(setdiff(campos, names(cl)), collapse = ', '))
}
if (!('id' %in% names(coh))) stop('El TSV de cohorte debe contener la columna id.')
if (!all(cl[['project.project_id']] == 'TCGA-BRCA')) {
  stop('ERROR: clinical.tsv contiene proyectos que NO son TCGA-BRCA.')
}

ids_coh <- unique(coh[['id']])
ids_cl <- unique(cl[['cases.case_id']])
if (anyNA(ids_coh) || anyDuplicated(coh[['id']])) {
  stop('La cohorte contiene IDs vacios o repetidos.')
}
if (!setequal(ids_coh, ids_cl)) {
  stop('Los Case ID de la cohorte y clinica no coinciden. ',
       'Cohorte sin clinica: ', sum(!(ids_coh %in% ids_cl)),
       ' | Clinica fuera de cohorte: ', sum(!(ids_cl %in% ids_coh)))
}

# IMPORTANTE: el TSV clinico tiene multiples filas por paciente debidas a
# diagnósticos, tratamientos y otras relaciones uno-a-muchos del GDC.
# Solo se aceptan estadios del diagnostico identificado como enfermedad primaria.
cl_primaria <- cl[tolower(as.character(
  cl[['diagnoses.diagnosis_is_primary_disease']])) == 'true', , drop = FALSE]

no_dato <- c('', "'--", '--', 'not reported', 'unknown',
             'not applicable', 'not evaluated')
limpiar_estadio <- function(v) {
  v <- trimws(as.character(v))
  v[is.na(v) | tolower(v) %in% no_dato] <- NA_character_
  v
}

# Una fila final por paciente. Si hay estadios incompatibles no elegimos uno
# arbitrariamente: se marcan como conflicto para revision.
resultado <- lapply(ids_coh, function(id) {
  sub <- cl_primaria[cl_primaria[['cases.case_id']] == id, , drop = FALSE]
  v <- unique(limpiar_estadio(sub[['diagnoses.ajcc_pathologic_stage']]))
  v <- v[!is.na(v)]
  utiles <- setdiff(v, 'Stage X')
  estado <- if (!nrow(sub)) {
    'Sin diagnostico primario'
  } else if (length(utiles) > 1L) {
    'Conflicto entre estadios primarios'
  } else if (length(utiles) == 1L) {
    'Disponible'
  } else if ('Stage X' %in% v) {
    'Stage X'
  } else {
    'Sin dato'
  }
  estadio <- if (estado == 'Disponible') utiles[[1]] else NA_character_
  data.frame(case_id = id, estadio_ajcc = estadio,
             estado_estadio = estado, stringsAsFactors = FALSE)
})
pacientes <- do.call(rbind, resultado)
stopifnot(nrow(pacientes) == length(ids_coh), !anyDuplicated(pacientes$case_id))

# Una primera informacion demografica por paciente, sin inventar imputaciones.
# Se conserva solamente en data/processed/, excluido por .gitignore.
cl_una <- cl[!duplicated(cl[['cases.case_id']]), , drop = FALSE]
ix <- match(pacientes$case_id, cl_una[['cases.case_id']])
pacientes$paciente_tcga <- cl_una[['cases.submitter_id']][ix]
raw_edad <- as.character(cl_una[['demographic.age_at_index']][ix])
raw_edad[tolower(trimws(raw_edad)) %in% no_dato] <- NA_character_
pacientes$edad_anos <- suppressWarnings(as.numeric(raw_edad))
if ('demographic.gender' %in% names(cl_una)) {
  pacientes$genero <- cl_una[['demographic.gender']][ix]
}

# Agrupacion AJCC para informe. No se asigna Stage X a IV ni se mezclan
# estadios patologicos con estadios clinicos.
agrupar <- function(s) {
  if (is.na(s)) return('Sin estadio utilizable')
  if (s %in% c('Stage I', 'Stage IA', 'Stage IB', 'Stage IC')) return('I')
  if (s %in% c('Stage II', 'Stage IIA', 'Stage IIB', 'Stage IIC')) return('II')
  if (s %in% c('Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC')) return('III')
  if (s == 'Stage IV') return('IV')
  if (s %in% c('Stage 0', 'Stage 0is')) return('0')
  'Otro estadio: revisar'
}
pacientes$estadio_agrupado <- vapply(pacientes$estadio_ajcc,
                                      agrupar, character(1))

resumen <- data.frame(
  indicador = c('Filas de clinical.tsv',
                'Pacientes clinicos unicos',
                'Pacientes en cohorte GDC',
                'Pacientes con diagnostico primario',
                'Pacientes con estadio patologico utilizable',
                'Pacientes sin estadio patologico utilizable',
                'Pacientes con estadios primarios contradictorios',
                'Pacientes con edad disponible'),
  n = c(nrow(cl), length(ids_cl), nrow(pacientes),
        sum(pacientes$estado_estadio != 'Sin diagnostico primario'),
        sum(pacientes$estado_estadio == 'Disponible'),
        sum(pacientes$estado_estadio != 'Disponible'),
        sum(pacientes$estado_estadio == 'Conflicto entre estadios primarios'),
        sum(!is.na(pacientes$edad_anos))),
  stringsAsFactors = FALSE
)

niveles <- c('0', 'I', 'II', 'III', 'IV',
             'Sin estadio utilizable', 'Otro estadio: revisar')
frecuencias <- table(factor(pacientes$estadio_agrupado, levels = niveles))
distribucion <- data.frame(estadio = names(frecuencias),
                           pacientes = as.integer(frecuencias),
                           stringsAsFactors = FALSE)
distribucion <- distribucion[distribucion$pacientes > 0, , drop = FALSE]

estados <- as.data.frame(table(pacientes$estado_estadio),
                         stringsAsFactors = FALSE)
names(estados) <- c('estado_estadio', 'pacientes')

utils::write.csv(resumen, file.path(salidas, '01_resumen_clinico.csv'),
                 row.names = FALSE)
utils::write.csv(distribucion, file.path(salidas, '02_estadios_ajcc.csv'),
                 row.names = FALSE)
utils::write.csv(estados, file.path(salidas, '03_faltantes_estadios.csv'),
                 row.names = FALSE)
utils::write.csv(pacientes,
                 file.path(salida_local, 'clinica_por_paciente.csv'),
                 row.names = FALSE)

nota <- c(
  '# TCGA-BRCA: disponibilidad de datos clinicos', '',
  'Fuente: GDC, exportacion de TCGA-BRCA del 2026-11-16.',
  'Entradas locales: data/raw/TCGA_BRCA/clinical/clinical.tsv y',
  'data/raw/TCGA_BRCA/cohort_TCGA_BRCA_LSQ.2026-11-16.tsv.', '',
  '## Metodo reproducible',
  '- Verificacion estricta de proyecto TCGA-BRCA y concordancia de Case ID.',
  '- Multiples filas de clinical.tsv se consolidan por paciente.',
  '- Estadio AJCC obtenido SOLO de diagnosticos con diagnosis_is_primary_disease=true.',
  '- Stage X, sin dato y ausencia de diagnostico primario: no utilizables.',
  '- Estadios contradictorios del diagnostico primario se marcan como conflicto.',
  '- Se reportan estadios I-IV; la agrupacion conserva un apartado para 0 u otros.',
  '- Los datos de paciente por identificador se dejan en data/processed/ (no GitHub).', '',
  '## Alcance',
  'Este script valida METADATOS clinicos. No descarga matrices STAR Counts,',
  'no realiza deconvolucion ni afirma asociacion entre estadios y proporciones.',
  'Antes de analizar RNA-seq bulk debe descargarse la subcohorte seleccionada',
  'y comprobarse la correspondencia de genes con la referencia scRNA-seq.'
)
writeLines(nota, file.path(salidas, '04_metodologia_clinica.md'),
           useBytes = TRUE)

cat('\n=== TCGA-BRCA | Validacion clinica terminada ===\n')
print(resumen, row.names = FALSE)
cat('\nDistribucion de estadios AJCC:\n')
print(distribucion, row.names = FALSE)
cat('\nArchivos publicables: ', salidas, '\n', sep = '')
cat('Tabla paciente-nivel LOCAL: ', salida_local, '\n', sep = '')
