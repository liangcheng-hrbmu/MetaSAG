# MetaSAG Usage 
## Step 4. Droplet clustering of potentially unknown species

## Func 1：FilterUnassignedCells(AnnotationFile, inputFastqDir, OutputDir)
- **Function Description:**

Automatically selects droplets that potentially belong to unknown species based on upstream annotation files and generates independent working directories for their aggregation.

- **Required Parameters:**
```
AnnotationFile  --      The path to the `CellAnno.txt` file generated in the upstream `MetaSAG_MetaPhlAnAsign` step. 

inputFastq      --      The directory containing the individual FASTQ files for each cell.

outputDir       --      Path to the output directory.
```


## Func 2：ClusterSAG(inputFastq,outputDir)
- **Function Description:**

Performs clustering (first round) on selected potentially unknown species.

- **Required Parameters:**
```

inputFastq      --      Path to the fastq files directory of potentially unknown species droplets.

outputDir       --      Path to save the clustering results.

```

- **Optional Parameters:**
```
ReadsEnd        --      Type of droplet sequencing files (single-end or paired-end).
                        Default: single-end, ReadsEnd='Single'.
                        For paired-end, set ReadsEnd='Pair'.
                        
SpadesEnv       --      Conda environment required for running Spades.py.
                        Default: None
                        
SourmashEnv     --      Conda environment required for running Sourmash.
                        Default: None            

ClusterThreshold--      The cut-off threshold used for forming flat clusters.
                        Default: 0.95.  

ClusterCriterion--      The criterion used to determine cluster formation.
                        Available options are 'distance' and 'inconsistent'.
                        Default: 'distance'.     
```



## Func 3：ClusterBin(OldRoundFold,NewRoundFold)

- **Function Description:**

Performs an additional round of clustering based on previous clustering results as needed.


- **Required Parameters:**
```

OldRoundFold    --      Root path of the previous clustering results (input data for the function).

outputDir       --      Path to save the results of the additional clustering.

```

- **Optional Parameters:**
```
ReadsEnd        --      Type of droplet sequencing files (single-end or paired-end).
                        Default: single-end, ReadsEnd='Single'.
                        For paired-end, set ReadsEnd='Pair'.
                        
SpadesEnv       --      Conda environment required for running Spades.py.
                        Default: None
                        
SourmashEnv     --      Conda environment required for running Sourmash.
                        Default: None

ClusterThreshold--      The cut-off threshold used for forming flat clusters.
                        Default: 0.95.  

ClusterCriterion--      The criterion used to determine cluster formation.
                        Available options are 'distance' and 'inconsistent'.
                        Default: 'distance'.                      
```



```bash

# Execution Command Examples

from MetaSAG import UnknownSAG as usag

fastqDir = Target_Path + 'CellBarn/'  

cellAnno =  Target_Path + 'MetaPhlAnAsign/MPAsign/CellAnno.txt'

outputdir = Target_Path + 'UnknownSAG/'

fastq_outputdir = outputdir + "Fastq/"

usag.FilterUnassignedCells(cellAnno, fastqDir, fastq_outputdir)

usag.ClusterSAG(fastq_outputdir,outputdir,ReadsEnd = 'Pair',SourmashEnv='sourmash') # Initial clustering of unknown species droplets

Round1Dir=os.path.join(outputdir,'Round1')

Round2Dir=os.path.join(outputdir,'Round2')

usag.ClusterBin(Round1Dir,Round2Dir,ReadsEnd = 'Pair',SourmashEnv='sourmash') # Re-clustering of unknown species droplet bins based on initial results

```

## Note on Test Data Limitation

This step is designed for droplets that remain unassigned after the MetaPhlAn4-based cell annotation step.

It uses the `CellAnno.txt` file generated in Step 3 to select cells labeled as `UnAsignedCell`, then clusters and assembles these droplets to explore potentially unknown species.

Because the Step 1–5 test dataset is intentionally small and mainly used to verify that the workflow can run, it may contain too few unassigned droplets, or the selected droplets may have insufficient sequencing depth. In that case, Step 4 may not generate meaningful clustering or assembly results.

For real biological interpretation, users should run this step on their own data with enough unassigned droplets and sufficient read depth.

Required input from previous steps:

```text
Target_Path/MetaPhlAnAsign/MPAsign/
└── CellAnno.txt

Target_Path/Barn/Cell/
├── Cell<SampleID><BarcodeIndex>_R1.fastq
└── Cell<SampleID><BarcodeIndex>_R2.fastq
```

`CellAnno.txt` is the key input file used to identify `UnAsignedCell` droplets. The corresponding cell FASTQ files are then selected from `Barn/Cell/` for unknown-species clustering.

## Usage Template

```python
from MetaSAG import UnknownSAG as usag
import os

Target_Path = "Your/Result/Path/"

cellAnno = Target_Path + "MetaPhlAnAsign/MPAsign/CellAnno.txt"
fastqDir = Target_Path + "Barn/Cell/"
outputDir = Target_Path + "UnknownSAG/"
fastqOutputDir = outputDir + "Fastq/"

usag.FilterUnassignedCells(cellAnno, fastqDir, fastqOutputDir)

usag.ClusterSAG(
    fastqOutputDir,
    outputDir,
    ReadsEnd="Pair",
    SpadesEnv=SpadesEnv,
    SourmashEnv="sourmash"
)

Round1Dir = os.path.join(outputDir, "Round1")
Round2Dir = os.path.join(outputDir, "Round2")

usag.ClusterBin(
    Round1Dir,
    Round2Dir,
    ReadsEnd="Pair",
    SpadesEnv=SpadesEnv,
    SourmashEnv="sourmash"
)
```

## Expected Output

If enough unassigned droplets are available, this step may generate:

```text
Target_Path/UnknownSAG/
├── Fastq/                         # FASTQ files of selected UnAsignedCell droplets
├── Round1/
│   ├── CellFasta/                 # SPAdes assemblies for selected cells
│   ├── sig/                       # sourmash signatures
│   └── round1_cluster_cell.txt    # Initial clustering result
└── Round2/                        # Optional re-clustering result
```

For small test datasets, these outputs may be empty or not biologically meaningful.