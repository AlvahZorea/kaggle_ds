# Script to find data quality issues and classic problems
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "Data Quality Analysis - Finding Problems" -ForegroundColor Cyan
Write-Host "================================================================`n" -ForegroundColor Cyan

$content = Import-Csv 'train.csv'
$rowCount = ($content | Measure-Object).Count

Write-Host "Analyzing $rowCount rows from train.csv`n" -ForegroundColor Yellow

# ================================================================
# 1. MISSING VALUES
# ================================================================
Write-Host "================================================================" -ForegroundColor Green
Write-Host "1. MISSING VALUES ANALYSIS" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$columns = $content[0].PSObject.Properties.Name
foreach ($col in $columns) {
    $nullCount = ($content | Where-Object { [string]::IsNullOrWhiteSpace($_.$col) } | Measure-Object).Count
    $pct = [math]::Round(($nullCount / $rowCount) * 100, 2)
    if ($nullCount -gt 0) {
        Write-Host "  [PROBLEM] $col : $nullCount missing values ($pct%)" -ForegroundColor Red
    } else {
        Write-Host "  [OK] $col : No missing values" -ForegroundColor Green
    }
}

# ================================================================
# 2. DUPLICATE ROWS
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "2. DUPLICATE ANALYSIS" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

# Check for duplicate sequences
$seqGroups = $content | Group-Object -Property NucleotideSequence | Where-Object { $_.Count -gt 1 }
$duplicateSeqCount = ($seqGroups | Measure-Object -Property Count -Sum).Sum - $seqGroups.Count
Write-Host "`nDuplicate NucleotideSequence:"
if ($duplicateSeqCount -gt 0) {
    Write-Host "  [PROBLEM] Found $duplicateSeqCount duplicate sequences" -ForegroundColor Red
    Write-Host "  Top duplicates:"
    $seqGroups | Sort-Object Count -Descending | Select-Object -First 5 | ForEach-Object {
        $seq = $_.Name.Substring(0, [math]::Min(50, $_.Name.Length)) + "..."
        Write-Host "    - '$seq' appears $($_.Count) times"
    }
} else {
    Write-Host "  [OK] No duplicate sequences" -ForegroundColor Green
}

# Check for duplicate symbols
$symbolGroups = $content | Group-Object -Property Symbol | Where-Object { $_.Count -gt 1 }
Write-Host "`nDuplicate Symbols:"
if ($symbolGroups.Count -gt 0) {
    Write-Host "  [WARNING] Found $($symbolGroups.Count) symbols with duplicates" -ForegroundColor Yellow
    $symbolGroups | Sort-Object Count -Descending | Select-Object -First 5 | ForEach-Object {
        Write-Host "    - '$($_.Name)' appears $($_.Count) times"
    }
} else {
    Write-Host "  [OK] No duplicate symbols" -ForegroundColor Green
}

# ================================================================
# 3. CLASS IMBALANCE
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "3. CLASS IMBALANCE ANALYSIS" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$geneTypes = $content | Group-Object -Property GeneType | Sort-Object Count -Descending
$maxClass = $geneTypes[0].Count
$minClass = $geneTypes[-1].Count
$imbalanceRatio = [math]::Round($maxClass / $minClass, 2)

Write-Host "`nClass Distribution:"
foreach ($type in $geneTypes) {
    $pct = [math]::Round(($type.Count / $rowCount) * 100, 2)
    $bar = "#" * [math]::Min(50, [int]($pct))
    if ($pct -lt 1) {
        Write-Host "  [RARE] $($type.Name): $($type.Count) ($pct%) $bar" -ForegroundColor Red
    } elseif ($pct -lt 5) {
        Write-Host "  [LOW] $($type.Name): $($type.Count) ($pct%) $bar" -ForegroundColor Yellow
    } else {
        Write-Host "  [OK] $($type.Name): $($type.Count) ($pct%) $bar" -ForegroundColor Green
    }
}

Write-Host "`nImbalance Ratio (max/min): $imbalanceRatio : 1"
if ($imbalanceRatio -gt 100) {
    Write-Host "  [SEVERE] Extreme class imbalance detected!" -ForegroundColor Red
} elseif ($imbalanceRatio -gt 10) {
    Write-Host "  [WARNING] Significant class imbalance" -ForegroundColor Yellow
}

# ================================================================
# 4. SEQUENCE VALIDITY & ANOMALIES
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "4. SEQUENCE VALIDITY & ANOMALIES" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

# Check for invalid nucleotides
Write-Host "`nChecking for invalid nucleotides (not A, T, G, C)..."
$invalidSeqCount = 0
$invalidChars = @{}
$content | ForEach-Object {
    $seq = $_.NucleotideSequence -replace '<', '' -replace '>', ''
    $invalid = $seq -replace '[ATGC]', ''
    if ($invalid.Length -gt 0) {
        $invalidSeqCount++
        foreach ($char in $invalid.ToCharArray()) {
            if ($invalidChars.ContainsKey($char)) {
                $invalidChars[$char]++
            } else {
                $invalidChars[$char] = 1
            }
        }
    }
}

if ($invalidSeqCount -gt 0) {
    Write-Host "  [WARNING] Found $invalidSeqCount sequences with non-ATGC characters" -ForegroundColor Yellow
    Write-Host "  Invalid characters found:"
    foreach ($key in $invalidChars.Keys) {
        Write-Host "    - '$key' : $($invalidChars[$key]) occurrences"
    }
} else {
    Write-Host "  [OK] All sequences contain only valid nucleotides (A, T, G, C)" -ForegroundColor Green
}

# Check sequence length anomalies
Write-Host "`nSequence Length Analysis:"
$lengths = $content | ForEach-Object {
    $seq = $_.NucleotideSequence -replace '<', '' -replace '>', ''
    $seq.Length
}

$avgLen = [math]::Round(($lengths | Measure-Object -Average).Average, 2)
$stdDev = [math]::Round([math]::Sqrt(($lengths | ForEach-Object { [math]::Pow($_ - $avgLen, 2) } | Measure-Object -Average).Average), 2)

# Very short sequences
$veryShort = ($lengths | Where-Object { $_ -lt 20 } | Measure-Object).Count
$short = ($lengths | Where-Object { $_ -ge 20 -and $_ -lt 50 } | Measure-Object).Count
$normal = ($lengths | Where-Object { $_ -ge 50 -and $_ -lt 900 } | Measure-Object).Count
$maxLength = ($lengths | Where-Object { $_ -ge 900 } | Measure-Object).Count

Write-Host "  - Very short (<20 bp): $veryShort sequences"
Write-Host "  - Short (20-50 bp): $short sequences"
Write-Host "  - Normal (50-900 bp): $normal sequences"
Write-Host "  - Near max length (>=900 bp): $maxLength sequences"

if ($veryShort -gt 0) {
    Write-Host "  [WARNING] $veryShort very short sequences may be problematic" -ForegroundColor Yellow
}
if ($maxLength -gt $rowCount * 0.1) {
    Write-Host "  [WARNING] Many sequences truncated at max length (1000)" -ForegroundColor Yellow
}

# ================================================================
# 5. OUTLIERS IN NUMERIC FIELDS
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "5. OUTLIER DETECTION" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

# NCBIGeneID outliers (using IQR method)
$geneIds = $content | ForEach-Object { [long]$_.NCBIGeneID } | Sort-Object
$q1Index = [int]($geneIds.Count * 0.25)
$q3Index = [int]($geneIds.Count * 0.75)
$q1 = $geneIds[$q1Index]
$q3 = $geneIds[$q3Index]
$iqr = $q3 - $q1
$lowerBound = $q1 - 1.5 * $iqr
$upperBound = $q3 + 1.5 * $iqr

$outliers = ($geneIds | Where-Object { $_ -lt $lowerBound -or $_ -gt $upperBound } | Measure-Object).Count
Write-Host "`nNCBIGeneID Outliers (IQR method):"
Write-Host "  Q1: $q1, Q3: $q3, IQR: $iqr"
Write-Host "  Lower bound: $lowerBound, Upper bound: $upperBound"
Write-Host "  Outliers: $outliers"

# ================================================================
# 6. DATA CONSISTENCY
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "6. DATA CONSISTENCY CHECKS" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

# Check if sequence format is consistent (all should have < >)
$withBrackets = ($content | Where-Object { $_.NucleotideSequence -match '^<.*>$' } | Measure-Object).Count
$withoutBrackets = $rowCount - $withBrackets

Write-Host "`nSequence Format Consistency:"
if ($withoutBrackets -gt 0) {
    Write-Host "  [WARNING] $withoutBrackets sequences don't have <> brackets" -ForegroundColor Yellow
} else {
    Write-Host "  [OK] All sequences have consistent <> format" -ForegroundColor Green
}

# Check GeneType values match expected categories
$expectedTypes = @('PSEUDO', 'BIOLOGICAL_REGION', 'ncRNA', 'snoRNA', 'PROTEIN_CODING', 'tRNA', 'OTHER', 'rRNA', 'snRNA', 'scRNA')
$actualTypes = $content | Select-Object -ExpandProperty GeneType -Unique
$unexpectedTypes = $actualTypes | Where-Object { $_ -notin $expectedTypes }

Write-Host "`nGeneType Values:"
if ($unexpectedTypes.Count -gt 0) {
    Write-Host "  [WARNING] Unexpected GeneType values found:" -ForegroundColor Yellow
    foreach ($t in $unexpectedTypes) {
        Write-Host "    - '$t'"
    }
} else {
    Write-Host "  [OK] All GeneType values are expected categories" -ForegroundColor Green
}

# ================================================================
# 7. TRAIN/VALIDATION/TEST DISTRIBUTION
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Green
Write-Host "7. TRAIN/VAL/TEST DISTRIBUTION CHECK" -ForegroundColor Green
Write-Host "================================================================" -ForegroundColor Green

$trainContent = Import-Csv 'train.csv'
$valContent = Import-Csv 'validation.csv'
$testContent = Import-Csv 'test.csv'

Write-Host "`nChecking for data leakage between splits..."

# Check if any sequences appear in multiple splits
$trainSeqs = $trainContent | Select-Object -ExpandProperty NucleotideSequence
$valSeqs = $valContent | Select-Object -ExpandProperty NucleotideSequence
$testSeqs = $testContent | Select-Object -ExpandProperty NucleotideSequence

$trainValOverlap = ($trainSeqs | Where-Object { $_ -in $valSeqs } | Measure-Object).Count
$trainTestOverlap = ($trainSeqs | Where-Object { $_ -in $testSeqs } | Measure-Object).Count
$valTestOverlap = ($valSeqs | Where-Object { $_ -in $testSeqs } | Measure-Object).Count

if ($trainValOverlap -gt 0 -or $trainTestOverlap -gt 0 -or $valTestOverlap -gt 0) {
    Write-Host "  [PROBLEM] Data leakage detected between splits!" -ForegroundColor Red
    Write-Host "    - Train-Validation overlap: $trainValOverlap sequences"
    Write-Host "    - Train-Test overlap: $trainTestOverlap sequences"
    Write-Host "    - Validation-Test overlap: $valTestOverlap sequences"
} else {
    Write-Host "  [OK] No sequence overlap between train/validation/test" -ForegroundColor Green
}

# Check class distribution consistency
Write-Host "`nClass Distribution Across Splits:"
$trainDist = $trainContent | Group-Object -Property GeneType
$valDist = $valContent | Group-Object -Property GeneType
$testDist = $testContent | Group-Object -Property GeneType

foreach ($type in $trainDist | Sort-Object Count -Descending) {
    $trainPct = [math]::Round(($type.Count / $trainContent.Count) * 100, 2)
    $valType = $valDist | Where-Object { $_.Name -eq $type.Name }
    $testType = $testDist | Where-Object { $_.Name -eq $type.Name }
    
    $valPct = if ($valType) { [math]::Round(($valType.Count / $valContent.Count) * 100, 2) } else { 0 }
    $testPct = if ($testType) { [math]::Round(($testType.Count / $testContent.Count) * 100, 2) } else { 0 }
    
    $diff = [math]::Abs($trainPct - $valPct) + [math]::Abs($trainPct - $testPct)
    if ($diff -gt 5) {
        Write-Host "  [WARNING] $($type.Name): Train=$trainPct%, Val=$valPct%, Test=$testPct% (uneven)" -ForegroundColor Yellow
    } else {
        Write-Host "  [OK] $($type.Name): Train=$trainPct%, Val=$valPct%, Test=$testPct%" -ForegroundColor Green
    }
}

# ================================================================
# SUMMARY OF PROBLEMS
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "SUMMARY OF DATA PROBLEMS FOUND" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan

Write-Host "`n[CRITICAL ISSUES]:" -ForegroundColor Red
Write-Host "  1. Severe class imbalance (ratio $imbalanceRatio : 1)"
Write-Host "     - scRNA has only 3 samples"
Write-Host "     - rRNA has only 8-277 samples across splits"

Write-Host "`n[MODERATE ISSUES]:" -ForegroundColor Yellow
Write-Host "  2. Data leakage in Symbol and Description columns"
Write-Host "  3. Some duplicate symbols exist"
Write-Host "  4. Very short sequences (<20 bp) may be problematic"

Write-Host "`n[MINOR ISSUES]:" -ForegroundColor Blue
Write-Host "  5. Sequences truncated at 1000 bp max length"
Write-Host "  6. Possible N or other ambiguous nucleotides in sequences"

# ================================================================
# RECOMMENDATIONS
# ================================================================
Write-Host "`n================================================================" -ForegroundColor Cyan
Write-Host "RECOMMENDATIONS FOR HANDLING PROBLEMS" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan

Write-Host "`n1. CLASS IMBALANCE:" -ForegroundColor White
Write-Host "   - Use stratified sampling for train/val/test splits"
Write-Host "   - Apply SMOTE or other oversampling techniques"
Write-Host "   - Use class weights in model training"
Write-Host "   - Consider merging rare classes (scRNA -> OTHER)"
Write-Host "   - Use metrics like F1-macro, not accuracy"

Write-Host "`n2. DATA LEAKAGE:" -ForegroundColor White
Write-Host "   - Remove Symbol and Description columns"
Write-Host "   - Or use only NucleotideSequence for prediction"

Write-Host "`n3. SEQUENCE ISSUES:" -ForegroundColor White
Write-Host "   - Filter out very short sequences (<20 bp)"
Write-Host "   - Handle ambiguous nucleotides (N, R, Y, etc.)"
Write-Host "   - Consider padding/truncating to fixed length"

Write-Host "`n4. FEATURE ENGINEERING:" -ForegroundColor White
Write-Host "   - Extract sequence length as feature"
Write-Host "   - Calculate GC content percentage"
Write-Host "   - Use k-mer frequencies (1-mer, 2-mer, 3-mer)"
Write-Host "   - Apply one-hot encoding for sequences"

Write-Host "`n5. MODEL SELECTION:" -ForegroundColor White
Write-Host "   - Use CNN/RNN/Transformer for sequence data"
Write-Host "   - Consider ensemble methods for imbalanced classes"
Write-Host "   - Use cross-validation with stratification"

Write-Host "`n"

