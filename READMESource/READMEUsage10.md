# MetaSAG Usage 
## Step 10. SGB Strain Evolution Analysis Function

A simple R script for analyzing bacterial strain evolution from multiple SGB folders.

## What It Does

This function processes **SGB*_Result** folders and:

    -Extracts strain genotypes from SNPs

    -Builds phylogenetic trees

    -Estimates divergence times

    -Generates reports and plots


## Required Input

<mark>Put your SGB folders in this structure:</mark>

your_data_folder/
<p>├── SGB1234_Result/<br>
│   ├── StrainCells.txt       # Strain assignments<br>
│   └── SGB1234_SNPpd.txt     # SNP matrix<br>
├── SGB5678_Result/<br>
│   ├── StrainCells.txt<br>
│   └── SGB5678_SNPpd.txt<br>
└── ...<p>


## Quick Start

```bash
# Basic usage - just point to your SGB folders
strain_analysis_multi_sgb.R -i /path/to/your/SGB_folders

# With custom output directory
strain_analysis_multi_sgb.R -i /path/to/SGB_folders -o /my/output
```

## Common Options
```bash

# Skip already analyzed SGBs
strain_analysis_multi_sgb.R -i /data/SGBs --skip-existing

# Change genome size (default: 5,000,000 bp)
strain_analysis_multi_sgb.R -i /data/SGBs --genome-size 3000000

# Use custom mutation rates
strain_analysis_multi_sgb.R -i /data/SGBs \
  --mutation-rate-medium 3e-6 \
  --generation-time-medium 3

# Run with debug output
strain_analysis_multi_sgb.R -i /data/SGBs --debug
```

## Output

Results go to Strain_Analysis_Results/ (or your specified output folder):

<p>Strain_Analysis_Results/<br>
├── SGB1234/<br>
│   ├── basic_statistics.tsv          # Basic counts<br>
│   ├── phylogeny_tree.png            # Strain tree<br>
│   ├── divergence_time.tsv           # Time estimates<br>
│   └── SGB1234_analysis_report.md    # Summary report<br>
├── SGB5678/<br>
└── batch_summary.tsv                 # All SGBs summary<p>

## Need Help?

Check the log files in the output folder, or run with --debug for more details.

Minimal command to get started:
```bash

strain_analysis_multi_sgb.R -i /your/data/folder
```
## Test Data

- [`Step_10_TestData`](../Example_data/Step_10_TestData)

The Step 1–5 small test dataset is mainly used to verify that the workflow can run and is not sufficient to support reliable strain evolution analysis. Therefore, Step 10 uses this independent test dataset to test SGB strain evolution analysis.

The test data include multiple SGB result folders. Each folder contains strain assignment information and an SNP matrix for one SGB.

```text
Step_10_TestData/
├── SGB1814_Result/
│   ├── StrainCells.txt
│   └── SGB1814_SNPpd.txt
├── SGB4269_Result/
│   ├── StrainCells.txt
│   └── SGB4269_SNPpd.txt
├── SGB4563_Result/
│   ├── StrainCells.txt
│   └── SGB4563_SNPpd.txt
└── ...
```

## Test Usage

```bash
strain_analysis_multi_sgb.R \
  -i ../Example_data/Step_10_TestData \
  -o Your/Result/Path/Strain_Analysis_Results
```

To skip SGB folders that have already been analyzed:

```bash
strain_analysis_multi_sgb.R \
  -i ../Example_data/Step_10_TestData \
  -o Your/Result/Path/Strain_Analysis_Results \
  --skip-existing
```

## Expected Output

```text
Your/Result/Path/Strain_Analysis_Results/
├── SGB1814/
│   ├── basic_statistics.tsv
│   ├── phylogeny_tree.png
│   ├── divergence_time.tsv
│   └── SGB1814_analysis_report.md
├── SGB4269/
├── SGB4563/
└── batch_summary.tsv
```

The exact set of SGB output folders depends on the valid `SGB*_Result` folders found in the input directory.
