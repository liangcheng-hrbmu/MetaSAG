# MetaSAG Usage 
## Step 8. HUMAnN Path.

## Class：HP(FastqDir,ResultDir)
- **Class Function:**

Classifies all sequencing reads and annotates biological pathways using HUMAnN.

- **Required Parameters:**
```

FastqDir        --      Path to the fastq files directory.

ResultDir       --      Path to save the results.

```

## Func 1：Diamond()

- **Function Description:**

Perform Diamond alignment on the FASTQ files under the FastqDir directory.


- **Optional Parameters:**
```
DiamondDB       --      Path to the Diamond reference database.
                        Default: None
                        
Diamondenv      --      Conda environment required for running Diamond.
                        Default: None

```



## Func 2：Uniref2Matrix()

- **Function Description:**

Generates a read count matrix of corresponding Uniref segments in each cell based on Diamond alignment results.


- **Optional Parameters:**

```

MinUnirefNum        --      Minimum count threshold for annotated Uniref reads in a cell.
                            Default: 5
    
MinCellNum          --      Minimum count threshold for the number of cells where Uniref-mapped reads appear.
                            Default: 5

```



## Func 3：SeuratCluster()

- **Function Description:**

Inputs the Cell-Uniref count matrix from Uniref2Matrix to cluster cells using Seurat, generating cell clusters related to Uniref counts.


![Seurat_Pheatmap](Seurat_Pheatmap.png)
![Seurat_Umap](Seurat_Umap.png)



## Func 4：HUMAnNPath(CellAnno, Group)

- **Function Description:**

Combines cell grouping information from CellAnno and Group to annotate biological pathways for each cell group using HUMAnN based on their Uniref annotations.


- **Required Parameters:**

```
CellAnno        --      Path to the cell grouping information file.

Group           --      Name of the grouping column in CellAnno.

```


- **可选参数：**
```
HUMAnNenv       --      Conda environment required for running HUMAnN.
                        Default: None

```


Eg. CellAnno (SeuratResult/KnownSGBCell_ClusterCell.txt)

| Cluster |     Cell      |
|:-------:|:-------------:|
|    1    | Sam1025_10012 |
|    1    | Sam1025_10168 |
|    2    | Sam1025_10335 |
|    2    | Sam1025_10713 |
|   ...   |      ...      |


- **Result:**

![HUMAnNPath](HUMAnNPath.png)



```

# Execution Command Examples

from MetaSAG import HUMAnNPath as hp

# Create an HP object

fastqDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/HUMAnNPath/input/fastq' #292Mb

resultDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/HUMAnNPath/result'

obj=hp.HP(fastqDir,resultDir)



# Perform Diamond alignment on each fastq file

obj.Diamond(DiamondDB='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4/Func_Anno/Humann3/DB/uniref/uniref90_201901b_full.dmnd')

# obj.DiamondDir 

# '/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/HUMAnNPath/result/DiamondDir'
# If the user performs Diamond alignment independently, modify obj.DiamondDir to specify the directory containing Diamond alignment results for subsequent analysis.
# obj.DiamondDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/HUMAnNPath/input/Diamond'
obj.Uniref2Matrix()
# Uniref2Matrix took 4.7900 seconds to execute.
obj.SeuratCluster()
# SeuratCluster took 16.6086 seconds to execute.

cellAnno='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData/testData/HUMAnNPath/result/SeuratResult/KnownSGBCell_ClusterCell.txt'
obj.HUMAnNPath(cellAnno,'Cluster',HUMAnNenv='humann')
# HUMAnNPath took 581.9383 seconds to execute.



```