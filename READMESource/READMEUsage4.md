# MetaSAG Usage 
## Step 4. Quality control and integration annotation of assembled genome.

## Func 1：FastaQC1(FastaDir,FastaOut)

- **Function Description:**
Removes overly short contigs from all FASTA files in the specified directory.

- **Required Parameters:**
```
FastaDir        --      Path to the directory containing uncorrected FASTA files.
                        Files must end with .fasta.

FastaOut        --      Path to save the directory of corrected FASTA files.

```````
- **Optional Parameters:**
```
minlen          --      Minimum contig length threshold in FASTA files.
                        Default: 500 bp.

```


## Func 2：FastaQC2(FastaDir,FastaOut,CheckmFile)

- **Function Description:**
Classifies corrected genome FASTA files into "qualified," "needs re-correction," or "abandon" based on CheckM results.

- **Required Parameters:**
```
FastaDir        --      Path to the directory containing corrected FASTA files.

FastaOut        --      Path to save the classification results of FASTA files.

CheckmFile      --      CheckM result table for the corrected FASTA files.

```````
- **Optional Parameters:**
```
FastgDir        --      If a fastg folder is provided, copy corresponding fastg files to the Bandage folder; 
                        if not, only place corresponding FASTA files in the Bandage folder.
                        Default: None.

```


## Func 3：Summary(FastaSGB,CheckmFile,GTDBFile)

- **Function Description:**

Integrates results from MetaPhlAnN4 classification files (FastaSGB), CheckM quality inspection files (CheckmFile), and GTDB-TK annotation files (GTDBFile) for assembled genomes.

- **Required Parameters:**
```
FastaSGB        --      Path to the MetaPhlAnN4 annotation classification file for assembled genomes.

CheckmFile      --      CheckM quality inspection file for assembled genomes.

GTDBFile        --      GTDB-TK annotation file for assembled genomes.

```


- **Optional Parameters:**
```
outputSummary   --      Path to save the integrated annotation results file.
                        Default: None (results are returned as a DataFrame object).

```



Eg. FastaSGB (first column: bin genome file name without .fasta suffix; second column: MetaPhlAn4 classification SGB ID)

|   Bin   |   SGB    |
|:-------:|:--------:|
| genome1 | SGB1024  |
| genome2 | SGB28348 |
| genome3 | SGB4348  |
|   ...   |   ...    |


- **Result:**

Eg. Summary.txt

|   Bin   |   SGB    |                                                         SGBName                                                          | Completeness | Contamination | closest_reference |                                                     closest_taxonomy                                                     |closest_ani|ContamFlag|CompleteFlag|QC|
|:-------:|:--------:|:------------------------------------------------------------------------------------------------------------------------:|:------------:|:-------------:|:-----------------:|:------------------------------------------------------------------------------------------------------------------------:|:---:|:---:|:---:|:---:|
| genome1 | SGB1024  |              k__Bacteria\|p__Bacteroidetes\|c__CFGB343\|o__OFGB343\|f__FGB343\|g__GGB781\|s__GGB781_SGB1024              |    96.02     |     13.66     |  GCA_000434395.1  |                d__Bacteria;p__Firmicutes;c__Bacilli;o__RFN20;f__CAG-288;g__CAG-568;s__CAG-568 sp000434395                |97.29|fail|pass|Low|
| genome2 | SGB28348 | k__Bacteria\|p__Bacteroidetes\|c__Cytophagia\|o__Cytophagales\|f__Cytophagaceae\|g__Aquirufa\|s__Aquirufa_antheringensis |    86.21     |     71.32     |        no         |                                                            no                                                            |no|fail|pass|Low|
| genome3 | SGB4348  |            k__Bacteria\|p__Firmicutes\|c__CFGB4806\|o__OFGB4806\|f__FGB4806\|g__GGB51647\|s__GGB51647_SGB4348            |    98.66     |     7.31      |  GCA_000432435.1  | d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Fimenecus;s__Fimenecus sp000432435 |98.92|pass|pass|Medium|
|   ...   |   ...    |                                                           ...                                                            |     ...      |      ...      |        ...        |                                                           ...                                                            |...|...|...|...|




```
#Execution Command Examples

from MetaSAG import BinQCAnno as bqa


# Remove contigs <500bp in length

FastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC1Fasta'

FastaOut='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/result/QC1Fasta'

bqa.FastaQC1(FastaDir,FastaOut)



# Evaluate genome assembly quality based on CheckM files and determine which assemblies need further correction

FastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/fasta'

FastgDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/fastg'

CheckmFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/qa_result'

FastaOut='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/result/QC2Fasta'

bqa.FastaQC2(FastaDir,FastaOut,CheckmFile,FastgDir=FastgDir) 



# Integrate MetaPhlAnN4 annotation classification, CheckM quality assessment, and GTDBTK annotation classification results

FastaSGB='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/FastaSGB'

CheckmFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/qa_result'

GTDBFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/gtdbtk.bac120.summary.tsv'

outputSummary='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/summary.txt'

summary=bqa.Summary(FastaSGB,CheckmFile,GTDBFile)




```

