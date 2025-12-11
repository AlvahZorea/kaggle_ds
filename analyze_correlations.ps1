# Script to analyze correlations and identify redundant variables
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "Correlation Analysis & Redundant Variable Detection" -ForegroundColor Cyan
Write-Host "================================================================`n" -ForegroundColor Cyan

# Load train data for analysis
$content = Import-Csv 'train.csv'
$rowCount = ($content | Measure-Object).Count

Write-Host "Analyzing $rowCount rows from train.csv`n" -ForegroundColor Yellow

# ================================================================
# 1. Check GeneGroupMethod - is it constant?
# ================================================================
Write-Host "================================================================" -ForegroundColor Green
Write-Host "1. GeneGroupMethod Analysis" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$geneGroupMethods = $content | Group-Object -Property GeneGroupMethod
Write-Host "`nUnique values in GeneGroupMethod: $($geneGroupMethods.Count)"
foreach ($method in $geneGroupMethods) {
    $pct = [math]::Round(($method.Count / $rowCount) * 100, 2)
    Write-Host "  - '$($method.Name)': $($method.Count) ($pct%)"
}

if ($geneGroupMethods.Count -eq 1) {
    Write-Host "`n>>> REDUNDANT: GeneGroupMethod is CONSTANT - provides NO information!" -ForegroundColor Red
}

# ================================================================
# 2. Check Index column (first unnamed column)
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "2. Index Column (H1) Analysis" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$indexValues = $content | Select-Object -ExpandProperty H1
$uniqueIndexes = ($indexValues | Select-Object -Unique | Measure-Object).Count
Write-Host "`nUnique index values: $uniqueIndexes out of $rowCount"
$isSequential = $true
for ($i = 0; $i -lt [math]::Min(100, $indexValues.Count); $i++) {
    if ([int]$indexValues[$i] -ne $i) {
        $isSequential = $false
        break
    }
}
if ($isSequential) {
    Write-Host "Index appears to be sequential (0, 1, 2, ...)"
    Write-Host "`n>>> REDUNDANT: Index column is just row number - provides NO predictive information!" -ForegroundColor Red
}

# ================================================================
# 3. Check NCBIGeneID uniqueness and correlation with other vars
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "3. NCBIGeneID Analysis" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$uniqueGeneIds = ($content | Select-Object -ExpandProperty NCBIGeneID -Unique | Measure-Object).Count
Write-Host "`nUnique NCBIGeneID values: $uniqueGeneIds out of $rowCount"
if ($uniqueGeneIds -eq $rowCount) {
    Write-Host "NCBIGeneID is 100% unique (acts as primary key)"
    Write-Host "`n>>> LIKELY REDUNDANT: NCBIGeneID is a unique identifier - not useful for prediction!" -ForegroundColor Yellow
}

# Check if NCBIGeneID ranges correlate with GeneType
Write-Host "`nChecking NCBIGeneID ranges by GeneType:"
$geneTypes = $content | Group-Object -Property GeneType
foreach ($type in $geneTypes | Sort-Object Count -Descending | Select-Object -First 5) {
    $ids = $type.Group | ForEach-Object { [long]$_.NCBIGeneID }
    $minId = ($ids | Measure-Object -Minimum).Minimum
    $maxId = ($ids | Measure-Object -Maximum).Maximum
    $avgId = [math]::Round(($ids | Measure-Object -Average).Average, 0)
    Write-Host "  $($type.Name): Min=$minId, Max=$maxId, Avg=$avgId"
}

# ================================================================
# 4. Check Symbol correlation with GeneType
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "4. Symbol Analysis & Correlation with GeneType" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

# Check symbol patterns by GeneType
Write-Host "`nSymbol patterns by GeneType:"
foreach ($type in $geneTypes | Sort-Object Count -Descending | Select-Object -First 5) {
    $symbols = $type.Group | Select-Object -ExpandProperty Symbol
    # Count symbols ending with 'P' (pseudogene indicator)
    $endingWithP = ($symbols | Where-Object { $_ -match 'P\d*$' } | Measure-Object).Count
    $pctP = [math]::Round(($endingWithP / $type.Count) * 100, 2)
    
    # Count LOC symbols
    $locSymbols = ($symbols | Where-Object { $_ -like 'LOC*' } | Measure-Object).Count
    $pctLoc = [math]::Round(($locSymbols / $type.Count) * 100, 2)
    
    # Count MIR symbols
    $mirSymbols = ($symbols | Where-Object { $_ -like 'MIR*' } | Measure-Object).Count
    $pctMir = [math]::Round(($mirSymbols / $type.Count) * 100, 2)
    
    Write-Host "  $($type.Name):"
    Write-Host "    - Ending with P (pseudogene): $endingWithP ($pctP%)"
    Write-Host "    - Starting with LOC: $locSymbols ($pctLoc%)"
    Write-Host "    - Starting with MIR: $mirSymbols ($pctMir%)"
}

Write-Host "`n>>> Symbol contains patterns that CORRELATE with GeneType!" -ForegroundColor Yellow
Write-Host ">>> This could be a 'data leakage' issue - Symbol encodes the target!" -ForegroundColor Yellow

# ================================================================
# 5. Check Description correlation with GeneType
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "5. Description Analysis & Correlation with GeneType" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

Write-Host "`nDescription patterns by GeneType:"
foreach ($type in $geneTypes | Sort-Object Count -Descending | Select-Object -First 5) {
    $descriptions = $type.Group | Select-Object -ExpandProperty Description
    
    # Check for 'pseudogene' in description
    $hasPseudo = ($descriptions | Where-Object { $_ -match 'pseudogene' } | Measure-Object).Count
    $pctPseudo = [math]::Round(($hasPseudo / $type.Count) * 100, 2)
    
    # Check for 'RNA' in description
    $hasRNA = ($descriptions | Where-Object { $_ -match 'RNA' } | Measure-Object).Count
    $pctRNA = [math]::Round(($hasRNA / $type.Count) * 100, 2)
    
    # Check for 'protein' in description
    $hasProtein = ($descriptions | Where-Object { $_ -match 'protein' } | Measure-Object).Count
    $pctProtein = [math]::Round(($hasProtein / $type.Count) * 100, 2)
    
    Write-Host "  $($type.Name):"
    Write-Host "    - Contains 'pseudogene': $hasPseudo ($pctPseudo%)"
    Write-Host "    - Contains 'RNA': $hasRNA ($pctRNA%)"
    Write-Host "    - Contains 'protein': $hasProtein ($pctProtein%)"
}

Write-Host "`n>>> Description contains keywords that DIRECTLY REVEAL GeneType!" -ForegroundColor Red
Write-Host ">>> This is SEVERE data leakage - Description essentially encodes the target!" -ForegroundColor Red

# ================================================================
# 6. Nucleotide Sequence Analysis
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "6. NucleotideSequence Analysis" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

Write-Host "`nSequence length by GeneType:"
foreach ($type in $geneTypes | Sort-Object Count -Descending | Select-Object -First 5) {
    $lengths = $type.Group | ForEach-Object {
        $seq = $_.NucleotideSequence -replace '<', '' -replace '>', ''
        $seq.Length
    }
    $avgLen = [math]::Round(($lengths | Measure-Object -Average).Average, 2)
    $minLen = ($lengths | Measure-Object -Minimum).Minimum
    $maxLen = ($lengths | Measure-Object -Maximum).Maximum
    Write-Host "  $($type.Name): Avg=$avgLen, Min=$minLen, Max=$maxLen"
}

# ================================================================
# SUMMARY
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "SUMMARY: Redundant Variables Identified" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan

Write-Host "`n[REDUNDANT - Remove]:" -ForegroundColor Red
Write-Host "  1. Index (H1) - Just row number, no predictive value"
Write-Host "  2. GeneGroupMethod - Constant value (100% 'NCBI Ortholog')"
Write-Host "  3. NCBIGeneID - Unique identifier, no predictive value"

Write-Host "`n[DATA LEAKAGE - Careful!]:" -ForegroundColor Yellow
Write-Host "  4. Symbol - Contains patterns that encode GeneType"
Write-Host "     (e.g., 'P' suffix for PSEUDO, 'MIR' for ncRNA)"
Write-Host "  5. Description - Contains explicit GeneType keywords"
Write-Host "     (e.g., 'pseudogene', 'RNA', 'protein coding')"

Write-Host "`n[USEFUL FOR PREDICTION]:" -ForegroundColor Green
Write-Host "  6. NucleotideSequence - The main feature for prediction"

Write-Host "`n[TARGET VARIABLE]:" -ForegroundColor Magenta
Write-Host "  7. GeneType - The label to predict"

Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "RECOMMENDATIONS" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan
Write-Host "`n1. Remove: Index, GeneGroupMethod, NCBIGeneID (redundant)"
Write-Host "2. Consider removing: Symbol, Description (data leakage)"
Write-Host "3. Use only: NucleotideSequence to predict GeneType"
Write-Host "4. If using Symbol/Description - be aware of potential overfitting!"
Write-Host "`n"

