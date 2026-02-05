

# SGB Strain Evolution Analysis Function

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
Rscript strain_analysis_multi_sgb.R -i /path/to/your/SGB_folders

# With custom output directory
Rscript strain_analysis_multi_sgb.R -i /path/to/SGB_folders -o /my/output
```

## Common Options
```bash

# Skip already analyzed SGBs
Rscript strain_analysis_multi_sgb.R -i /data/SGBs --skip-existing

# Change genome size (default: 5,000,000 bp)
Rscript strain_analysis_multi_sgb.R -i /data/SGBs --genome-size 3000000

# Use custom mutation rates
Rscript strain_analysis_multi_sgb.R -i /data/SGBs \
  --mutation-rate-medium 3e-6 \
  --generation-time-medium 3

# Run with debug output
Rscript strain_analysis_multi_sgb.R -i /data/SGBs --debug
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

Rscript strain_analysis_multi_sgb.R -i /your/data/folder
```
