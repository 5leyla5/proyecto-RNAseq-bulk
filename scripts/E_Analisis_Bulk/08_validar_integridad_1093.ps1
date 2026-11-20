$ErrorActionPreference = "Stop"

$raiz = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path

$manifest = Join-Path $raiz "docs\E_Analisis_Bulk\06_cohorte_completa\manifest_1093_TCGA_BRCA.txt"
$destino = Join-Path $raiz "data\raw\TCGA_BRCA\expresion_completa"
$salida = Join-Path $raiz "docs\E_Analisis_Bulk\06_cohorte_completa\integridad_1093.csv"

$datos = @(Import-Csv -LiteralPath $manifest -Delimiter "`t")

if ($datos.Count -ne 1093) {
    throw "El manifiesto no contiene 1093 archivos."
}

$contador = 0

$resultados = foreach ($f in $datos) {
    $contador++

    if ($contador % 50 -eq 0 -or $contador -eq 1) {
        Write-Host "Verificando $contador de $($datos.Count)..."
    }

    $ruta = Join-Path (Join-Path $destino $f.id) $f.filename

    $tamanoReal = $null
    $md5Real = ""
    $estado = "FALTANTE"

    if (Test-Path -LiteralPath $ruta -PathType Leaf) {
        $tamanoReal = (Get-Item -LiteralPath $ruta).Length

        if ($tamanoReal -ne [long]$f.size) {
            $estado = "TAMANO_INCORRECTO"
        }
        else {
            try {
                $md5Real = (Get-FileHash -LiteralPath $ruta -Algorithm MD5).Hash

                if ($md5Real -ieq $f.md5.Trim()) {
                    $estado = "MD5_OK"
                }
                else {
                    $estado = "MD5_INCORRECTO"
                }
            }
            catch {
                $estado = "ERROR_LECTURA"
            }
        }
    }

    [PSCustomObject]@{
        id = $f.id
        filename = $f.filename
        size_esperado = $f.size
        size_real = $tamanoReal
        md5_esperado = $f.md5
        md5_real = $md5Real
        estado = $estado
    }
}

$resultados | Export-Csv -LiteralPath $salida -NoTypeInformation -Encoding UTF8

Write-Host "`n=== RESULTADO DE INTEGRIDAD ==="

$resultados |
    Group-Object estado |
    Select-Object Name, Count |
    Format-Table -AutoSize

Write-Host "Reporte guardado en: $salida"

$correctos = @($resultados | Where-Object { $_.estado -eq "MD5_OK" }).Count

if ($correctos -ne 1093) {
    throw "Integridad pendiente: $correctos de 1093 archivos con MD5 correcto."
}

Write-Host "VALIDACION COMPLETA: 1093/1093 MD5_OK"
