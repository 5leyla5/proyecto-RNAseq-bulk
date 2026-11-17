# E/04a - Descarga reproducible de 50 matrices STAR Counts TCGA-BRCA.
# Ejecutar desde PowerShell: 
#   & "scripts/E_Analisis_Bulk/04_descargar_bulk_TCGA_BRCA.ps1"
# También admite una ruta de gdc-client.exe opcional:
#   & "scripts/E_Analisis_Bulk/04_descargar_bulk_TCGA_BRCA.ps1" -GdcExe "C:\ruta\gdc-client.exe"
# Este script DESCARGA, no sustituye la verificacion MD5 del script R E/04.
# GDC Data Transfer Tool usada en la ejecucion inicial: version 2.3.

param(
  [string]$GdcExe = ""
)

$ErrorActionPreference = 'Stop'
$raiz = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$manifest = Join-Path $raiz 'docs\E_Analisis_Bulk\03_piloto\03_manifest_piloto_50_FINAL.txt'
$destino = Join-Path $raiz 'data\raw\TCGA_BRCA\expresion_piloto'
$dirRegistro = Join-Path $raiz 'docs\E_Analisis_Bulk\04_descarga'
$registro = Join-Path $dirRegistro '03_procedencia_gdc_client.json'

if (-not (Test-Path -LiteralPath (Join-Path $raiz 'README.md'))) {
  throw "No se encontro la raiz esperada del repositorio."
}
if (-not (Test-Path -LiteralPath $manifest -PathType Leaf)) {
  throw "No existe el manifiesto final. Ejecuta antes E/03: $manifest"
}

$filas = @(Import-Csv -LiteralPath $manifest -Delimiter "`t")
if ($filas.Count -ne 50 -or (@($filas.id | Select-Object -Unique)).Count -ne 50) {
  throw "El manifiesto debe contener exactamente 50 IDs unicos. Detectados: $($filas.Count)."
}
if (@($filas | Where-Object { -not $_.filename -or -not $_.md5 -or -not $_.size }).Count -gt 0) {
  throw "El manifiesto carece de filename, md5 o size en al menos una fila."
}

# Primero buscar en PATH; si no esta, buscar en Descargas como se hizo en Windows.
if (-not [string]::IsNullOrWhiteSpace($GdcExe)) {
  $cliente = (Resolve-Path -LiteralPath $GdcExe -ErrorAction Stop).Path
} else {
  $comando = Get-Command 'gdc-client.exe' -ErrorAction SilentlyContinue
  if (-not $comando) { $comando = Get-Command 'gdc-client' -ErrorAction SilentlyContinue }
  if ($comando -and $comando.Source) {
    $cliente = $comando.Source
  } else {
    $descargas = Join-Path $HOME 'Downloads'
    $candidatos = @()
    if (Test-Path -LiteralPath $descargas) {
      $candidatos = @(Get-ChildItem -LiteralPath $descargas -Filter 'gdc-client.exe' -File -Recurse -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending)
    }
    if ($candidatos.Count -eq 0) {
      throw "No se encontro gdc-client.exe. Descarga la herramienta y usa -GdcExe 'RUTA_AL_EXE'."
    }
    $cliente = $candidatos[0].FullName
  }
}
if (-not (Test-Path -LiteralPath $cliente -PathType Leaf)) {
  throw "La ruta encontrada no es un ejecutable: $cliente"
}

$version = ((& $cliente --version) | Out-String).Trim()
if ($LASTEXITCODE -ne 0) { throw "gdc-client --version termino con error." }
if ($version -ne '2.3') {
  Write-Warning "La descarga original uso gdc-client 2.3; esta instalacion reporta $version. Se registrara la diferencia."
}

# Crucial: gdc-client no crea necesariamente la carpeta de destino antes de escribir.
New-Item -Path $destino -ItemType Directory -Force | Out-Null
New-Item -Path $dirRegistro -ItemType Directory -Force | Out-Null

$shaManifest = (Get-FileHash -LiteralPath $manifest -Algorithm SHA256).Hash
$shaCliente = (Get-FileHash -LiteralPath $cliente -Algorithm SHA256).Hash
$bytes = [long](($filas | Measure-Object -Property size -Sum).Sum)
$fechaInicio = (Get-Date).ToUniversalTime().ToString('o')

Write-Host "GDC Data Transfer Tool $version"
Write-Host "Archivos en manifiesto: $($filas.Count); bytes esperados: $bytes"
Write-Host 'Descargando solamente la subcohorte de 50 archivos...'

# Permite reanudar usando el MISMO manifiesto y destino si se interrumpe.
& $cliente download -m $manifest -d $destino
$codigoSalida = $LASTEXITCODE

# No publicar rutas personales en GitHub. El registro guarda hashes, version,
# rutas RELATIVAS y estado de la ejecucion, nunca la ruta absoluta del usuario.
$metadatos = [ordered]@{
  herramienta = 'GDC Data Transfer Tool'
  version_gdc_client = $version
  sha256_ejecutable = $shaCliente
  fecha_inicio_utc = $fechaInicio
  fecha_fin_utc = (Get-Date).ToUniversalTime().ToString('o')
  manifest_relativo = 'docs/E_Analisis_Bulk/03_piloto/03_manifest_piloto_50_FINAL.txt'
  sha256_manifest = $shaManifest
  numero_archivos = $filas.Count
  bytes_esperados = $bytes
  destino_relativo = 'data/raw/TCGA_BRCA/expresion_piloto'
  codigo_salida = $codigoSalida
  advertencia = 'La finalizacion de gdc-client NO demuestra integridad. Ejecutar E/04 R (50 MD5 y tamanos).'
}
($metadatos | ConvertTo-Json -Depth 3) | Set-Content -LiteralPath $registro -Encoding UTF8

if ($codigoSalida -ne 0) {
  throw "gdc-client termino con codigo $codigoSalida. Revisa el error y reanuda usando este mismo script."
}
Write-Host 'Transferencia finalizada. Ahora valida con R:'
Write-Host "source('scripts/E_Analisis_Bulk/04_verificar_descarga_TCGA_BRCA.R')"
Write-Host "Registro de procedencia: docs/E_Analisis_Bulk/04_descarga/03_procedencia_gdc_client.json"
