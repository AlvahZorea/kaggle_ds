# Script to create data analysis report
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "Data Analysis Report" -ForegroundColor Cyan
Write-Host "================================================================`n" -ForegroundColor Cyan

$files = @('train.csv', 'validation.csv', 'test.csv')
$results = @()

foreach ($file in $files) {
    if (Test-Path $file) {
        Write-Host "`n================================================================" -ForegroundColor Yellow
        Write-Host "Analyzing file: $file" -ForegroundColor Yellow
        Write-Host "================================================================`n" -ForegroundColor Yellow
        
        $content = Import-Csv $file
        $rowCount = ($content | Measure-Object).Count
        
        Write-Host "General Information:" -ForegroundColor Green
        Write-Host "  - Number of rows: $rowCount"
        Write-Host "  - Number of columns: $($content[0].PSObject.Properties.Count)"
        
        Write-Host "`nList of Variables (Columns):" -ForegroundColor Green
        $i = 1
        foreach ($prop in $content[0].PSObject.Properties.Name) {
            Write-Host "  $i. $prop"
            $i++
        }
        
        # Analyze GeneType
        if ($content[0].PSObject.Properties.Name -contains 'GeneType') {
            Write-Host "`nLabel Analysis (GeneType):" -ForegroundColor Green
            $geneTypes = $content | Group-Object -Property GeneType | Sort-Object Count -Descending
            Write-Host "  - Number of unique categories: $($geneTypes.Count)"
            Write-Host "`n  Label Distribution:"
            foreach ($type in $geneTypes) {
                $pct = [math]::Round(($type.Count / $rowCount) * 100, 2)
                Write-Host "    - $($type.Name): $($type.Count) ($pct%)"
            }
        }
        
        # Analyze GeneGroupMethod
        if ($content[0].PSObject.Properties.Name -contains 'GeneGroupMethod') {
            Write-Host "`nGeneGroupMethod Analysis:" -ForegroundColor Green
            $methods = $content | Group-Object -Property GeneGroupMethod | Sort-Object Count -Descending
            foreach ($method in $methods) {
                $pct = [math]::Round(($method.Count / $rowCount) * 100, 2)
                Write-Host "  - $($method.Name): $($method.Count) ($pct%)"
            }
        }
        
        # Analyze NCBIGeneID
        if ($content[0].PSObject.Properties.Name -contains 'NCBIGeneID') {
            Write-Host "`nNCBIGeneID Analysis:" -ForegroundColor Green
            $uniqueIds = ($content | Select-Object -ExpandProperty NCBIGeneID -Unique | Measure-Object).Count
            Write-Host "  - Number of unique IDs: $uniqueIds"
            Write-Host "  - Number of duplicates: $($rowCount - $uniqueIds)"
        }
        
        # Analyze Symbol
        if ($content[0].PSObject.Properties.Name -contains 'Symbol') {
            Write-Host "`nSymbol Analysis:" -ForegroundColor Green
            $uniqueSymbols = ($content | Select-Object -ExpandProperty Symbol -Unique | Measure-Object).Count
            Write-Host "  - Number of unique symbols: $uniqueSymbols"
            Write-Host "  - Number of duplicates: $($rowCount - $uniqueSymbols)"
        }
        
        # Analyze nucleotide sequence lengths
        if ($content[0].PSObject.Properties.Name -contains 'NucleotideSequence') {
            Write-Host "`nNucleotide Sequence Analysis:" -ForegroundColor Green
            $lengths = $content | ForEach-Object {
                $seq = $_.NucleotideSequence -replace '<', '' -replace '>', ''
                $seq.Length
            }
            $avgLength = [math]::Round(($lengths | Measure-Object -Average).Average, 2)
            $minLength = ($lengths | Measure-Object -Minimum).Minimum
            $maxLength = ($lengths | Measure-Object -Maximum).Maximum
            $sortedLengths = $lengths | Sort-Object
            $medianLength = [math]::Round($sortedLengths[[math]::Floor($sortedLengths.Count / 2)], 2)
            
            Write-Host "  - Average length: $avgLength nucleotides"
            Write-Host "  - Minimum length: $minLength nucleotides"
            Write-Host "  - Maximum length: $maxLength nucleotides"
            Write-Host "  - Median length: $medianLength nucleotides"
        }
        
        $results += [PSCustomObject]@{
            File = $file
            Rows = $rowCount
            Columns = $content[0].PSObject.Properties.Count
        }
    } else {
        Write-Host "`nWarning: File $file not found" -ForegroundColor Red
    }
}

# Compare files
if ($results.Count -gt 1) {
    Write-Host "`n================================================================" -ForegroundColor Cyan
    Write-Host "Comparison between data files" -ForegroundColor Cyan
    Write-Host "================================================================`n" -ForegroundColor Cyan
    $results | Format-Table -AutoSize
}

Write-Host "`nReport completed!" -ForegroundColor Green
