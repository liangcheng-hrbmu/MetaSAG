# MetaSAG Usage 
## Step 9. Droplet clustering of potentially unknown species

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



```

# Execution Command Examples

from MetaSAG import UnknownSAG as usag

fastqDir = Target_Path + 'CellBarn/'  #846M

cellAnno =  Target_Path + 'MetaPhlAnAsign/MPAsign/CellAnno.txt'

outputdir = Target_Path + 'UnknownSAG/'

fastq_outputdir = outputdir + "Fastq/"

usag.FilterUnassignedCells(cellAnno, fastqDir, fastq_outputdir)

usag.ClusterSAG(fastq_outputdir,outputdir,ReadsEnd = 'Pair',SourmashEnv='sourmash') # Initial clustering of unknown species droplets

Round1Dir=os.path.join(outputdir,'Round1')

Round2Dir=os.path.join(outputdir,'Round2')

usag.ClusterBin(Round1Dir,Round2Dir,ReadsEnd = 'Pair',SourmashEnv='sourmash') # Re-clustering of unknown species droplet bins based on initial results

```