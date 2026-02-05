# MetaSAG Usage 
## Step 9. Droplet clustering of potentially unknown species

## Func 1：ClusterSAG(inputFastq,outputDir)
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

```



## Func 2：ClusterBin(OldRoundFold,NewRoundFold)

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

```



```

# Execution Command Examples

from MetaSAG import UnknownSAG as usag

fastqDir = Target_Path + 'Unknown/fastq/'  #846M

resultDir = Target_Path + 'Unknown/result/'

usag.ClusterSAG(fastqDir,resultDir,SourmashEnv='sourmash') # Initial clustering of unknown species droplets

Round1Dir=os.path.join(resultDir,'Round1')

Round2Dir=os.path.join(resultDir,'Round2')

usag.ClusterBin(Round1Dir,Round2Dir) # Re-clustering of unknown species droplet bins based on initial results

```