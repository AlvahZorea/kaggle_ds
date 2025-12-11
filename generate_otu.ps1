# Generate Fake OTU Table with Relative Abundances
# 100 OTUs across 10 cows

$random = New-Object System.Random(42)

# Bacterial genera for realistic OTU names
$genera = @(
    'Bacteroides', 'Prevotella', 'Ruminococcus', 'Clostridium', 'Lactobacillus',
    'Bifidobacterium', 'Streptococcus', 'Enterococcus', 'Escherichia', 'Faecalibacterium',
    'Roseburia', 'Blautia', 'Coprococcus', 'Dorea', 'Lachnospira',
    'Eubacterium', 'Anaerostipes', 'Butyricicoccus', 'Oscillospira', 'Dialister',
    'Megasphaera', 'Veillonella', 'Selenomonas', 'Succinivibrio', 'Fibrobacter',
    'Treponema', 'Butyrivibrio', 'Pseudobutyrivibrio', 'Acetitomaculum', 'Succiniclasticum'
)

# Generate OTU names
$otuNames = @()
for ($i = 0; $i -lt 100; $i++) {
    $genus = $genera[$i % $genera.Count]
    $speciesNum = [Math]::Floor($i / $genera.Count) + 1
    $otuNames += "OTU_{0:D3}_{1}_sp{2}" -f ($i + 1), $genus, $speciesNum
}

# Generate sample names (cows)
$sampleNames = @()
for ($i = 1; $i -le 10; $i++) {
    $sampleNames += "Cow_{0:D2}" -f $i
}

# Function to generate Dirichlet-like distribution
function Get-DirichletSample {
    param([double[]]$alpha)
    
    $gamma = @()
    foreach ($a in $alpha) {
        # Gamma approximation using sum of exponentials
        $sum = 0.0
        for ($j = 0; $j -lt [Math]::Max(1, [Math]::Floor($a)); $j++) {
            $sum += -[Math]::Log(1 - $random.NextDouble())
        }
        $gamma += $sum + (-[Math]::Log(1 - $random.NextDouble()) * ($a - [Math]::Floor($a)))
    }
    
    $total = ($gamma | Measure-Object -Sum).Sum
    return $gamma | ForEach-Object { $_ / $total }
}

# Generate abundance data
$abundanceData = @{}

foreach ($sample in $sampleNames) {
    # Create alpha parameters (log-normal-like for realistic microbiome distribution)
    $alpha = @()
    for ($i = 0; $i -lt 100; $i++) {
        $baseAlpha = [Math]::Exp($random.NextDouble() * 3 - 2) * 0.5 + 0.01
        $sampleVariation = $random.NextDouble() * 0.5 + 0.5
        $alpha += $baseAlpha * $sampleVariation
    }
    
    # Generate sample
    $abundances = Get-DirichletSample -alpha $alpha
    $abundanceData[$sample] = $abundances
}

# Create CSV content
$csvContent = "OTU," + ($sampleNames -join ",") + "`n"

for ($i = 0; $i -lt 100; $i++) {
    $row = $otuNames[$i]
    foreach ($sample in $sampleNames) {
        $row += "," + [Math]::Round($abundanceData[$sample][$i], 8)
    }
    $csvContent += $row + "`n"
}

# Save to file
$csvContent | Out-File -FilePath "otu_table.csv" -Encoding UTF8 -NoNewline

Write-Host "OTU table saved to otu_table.csv" -ForegroundColor Green
Write-Host "Shape: 100 OTUs x 10 samples"

