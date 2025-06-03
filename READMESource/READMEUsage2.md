# Step 2. Filter low-quality cells.

## Class：BCFilter(inputfastq,outputdir)
- **Class Description:**

Droplet Filtration Based on Read Count Distribution in Sample.

- **Required Parameters:**
```
inputFastq  --  Location of short-read sequencing file(s).
                For single-end FASTQ: Filename must end with .fastq.
                For paired-end FASTQ: Filenames must end with _R1.fastq and _R2.fastq, provided as a list (e.g., ["sample_R1.fastq", "sample_R2.fastq"]).

outputdir   --  Specifies the file path where result files will be stored.
```````


## Func 1：CellCountStatistic()

- **Function Description:**
Calculates and reports the number of sequencing reads assigned to each cell barcode in the input FASTQ file(s).

- **Result:**

Eg. bcread.txt

|    cell    | read  |
|:----------:|:-----:|
| Cell500000 |  332  |
| Cell500001 |  131  |
| Cell500002 | 88281 |
| Cell500003 |  420  |
| Cell500004 | 1593  |
|    ...     |  ...  |




## Func 2：getMinReads()

- **Function Description:**
Automatically calculates the optimal minimum read count threshold for filtering out low-quality cells based on the read count distribution across all cells in a sample.

- **Optional Parameters:**
```
FirstDrv_xlim       --      Specifies the horizontal axis range for first derivative plots.
                            Default Value:[-1000,5000]
                            
SecondDrv_xlim      --      Specifies the horizontal axis range for second derivative plots.
                            Default Value:[-1000,5000]
                            
```


## Func 3：BCMinReads()

- **Function Description:**
Filters cellular barcodes (droplets) based on a minimum read count threshold, generating a validated barcode output file.




```
#Execution Command Examples

from MetaSAG import BCFilter as bcf

#S10_Tag.fastq 6.1G
obj=bcf.BCFilter('./testData/BCFilter/input/S10_Tag.fastq','./testData/BCFilter/result/')

obj.CellCountStatistic()
#CellCountStatistic took 21.4586 seconds to execute.

# obj.BC_Count:Stores per-cell read count statistics in a modifiable tabular format for quality control and downstream analysis.

obj.getMinReads()
#getMinReads took 3.0193 seconds to execute.

# obj.min_reads:Minimum Droplet Read Count Threshold (Configurable)

obj.BCMinReads()

#obj.BC_Count_Filter:Generates a curated list of cellular barcodes that meet or exceed the specified minimum read count threshold (obj.min_reads)


```


![1stDerivate](1stDerivate.png)
![2ndDerivate](2ndDerivate.png)
![CDF_BCReads](CDF_BCReads.png)





## Func 4：SamsBCFilterStack(Sam,SavedCell,AllCell,outputDir)

- **Function Description:**

- Performs comparative analysis of cellular yields across multiple samples before and after quality filtering, with automated generation of publication-ready stacked barplots.



- **Required Parameters:**
```
Sam       --  Sample Name List

AllCell   --  Per-Sample Cell Count Statistics

SaveCell  --  List of Retained Cell Counts per Sample After Cell Filtering

outputDir --  Specifies the directory path where all output files will be saved.

```

```
#Execution Command Examples

AllCell = [3132, 4557, 2321, 3433, 4500]

SavedCell = [2000,3000,1734,2222,2000]

Sam=['S5','S6','S7','S8','S10']

outputDir=’./testData/BCFilter/result/’

from MetaSAG import BCFilter as bcf

bcf.SamBCFilterStack(Sam,SavedCell,AllCell,'./testData/BCFilter/result/')

#SamsBCFilterStack took 0.1096 seconds to execute.

```

![StackCellsPicture](StackCellsPicture.png)
