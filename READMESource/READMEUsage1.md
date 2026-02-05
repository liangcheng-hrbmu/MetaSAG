# Step 1. Distribute the reads in the sample to a file of individual droplets.

## Func 1：SAGSplit(inputFastq,CellBarn)

- **Function Description:**

The SAGSplit() function processes input sequencing files (in FASTQ format) and distributes each read into the corresponding droplet file based on its barcode sequence. 

Barcode identification is performed by the built-in FindBarcode() function or a user-provided alternative.


- **Required Parameters:**
```
inputFastq      Location of short-read sequencing file(s).
                For single-end FASTQ: Filename must end with .fastq.
                For paired-end FASTQ: Filenames must end with _R1.fastq and _R2.fastq, provided as a list (e.g., ["sample_R1.fastq", "sample_R2.fastq"]).

CellBarn        Output directory for droplet-split files. Stores partitioned results after processing.
```

- **Optional Parameters:**
```
FindBarcode     Custom barcode identification function. Accepts a single read string as input.
                Returns either:
                • Matched barcode string (e.g., "ATCGGTCA")
                • Error code string for mismatches (e.g., "Err", "Len", "W1", "BC")

warning	list	Designated error codes from FindBarcode indicating:
                • Len: Barcode length mismatch
                • W1: Weak signal/quality
                • BC: Barcode contamination

filterWarning   Filter control:
                • When True, reads returning warning codes are excluded from droplet files
                • When False, all reads are processed regardless of warnings

                    
                    

```


- **Single-End Data Processing**
```
#Input File Examples
test.fastq:

@NB501288_516_HC2NTBGXB:1:11101:18376:1042#CGAGGCTG/1
AGAGGTNGGAGTGATTGCTTGTGACGCCTTTGCCTCACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGTCCTTAACCATCCTTGAATACCTCGCTTGCTATTTTTTGTGCTTCTTTCCTCAGATATTGTGCCGTCTCATACATGAATGGTCTTT
+
AAAAAE#EEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAEEEEEAEEEEEEEEEE6AEEEEAEEEAEEEEEE/EE/EEEEAEAAEE/<EEE<AE</EEE<AAEE6EEAAEEE<EAEEEEEAEAEAE/A/EEAA</6/<EEE/EEEE</</E6AA/
@NB501288_516_HC2NTBGXB:1:11101:20037:1043#CGAGGCTG/1
GTTTGTNTGAGTGATTGCTTGTGACGCCTTGTTTGTTTTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGATAAATACGTATAGTACGATCAAAAACGCAAGAATATATCCGATCGCCCCAAGCGCCGGAATGCCGAGTATTTTCGGACTCATATC
+
AAAAAE#EEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE6/EEEEEEAEEAAEEE/EEEEEEEEEEEEEEEEEEEEEEEEA<AEEEEAAEEE<A/EE<<<<<AE/E/A<A<A<EEA66AAEEEEEAE<<<E<A<A6AA/

```


```
#Execution Command Examples

from MetaSAG import BarcodeDeal as bcd

#inputFastq File size: 11MB
inputFastq = "/DATA/test.fastq" # Single-end fastq file for analysis 

#Target_Path: A user-defined directory for storing output results. The path must terminate with a trailing forward slash (e.g., /path/to/output/).
CellBarn = Target_Path + "Barn/CellBarn_single/" # Destination subdirectory for processed results (Ensure trailing "/" is present).

bcd.SAGSplit(inputFastq,CellBarn)


```

```
#Result

{'W1': 774, 'Len': 0, 'BC': 1533}
SAGSplit took 1.8001 seconds to execute.

```






```
#Output File Examples

#318.fastq
#Each read header in the file starts with @Barcode_ID:

@318:NB501288_516_HC2NTBGXB:1:11101:8220:5572#CGAGGCTG/1
GTTTGTTTGAGTGATTGCTTGTGACGCCTTCCTGACACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGCTTGTATACAATATGCTTATAGTATACTCATATTTTCCTTAAAAATCAATATTTTATCTCACGATTTTAAATCTGAATTTTCCATTT
+
AAAAA/EEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEEEEEEEEEEAEEEEEE6AEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE<AEAEAEEAE<<EEAEEEEE6AEEEE/A/E<<<///AEE66AEEEEAEAEEAE/A</<</A
@318:NB501288_516_HC2NTBGXB:1:11101:23509:16472#CGAGGCTG/1
GTTTGTTTGAGTGATTGCTTGTGACGCCTTCCTGACACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGTACAGCCGCATTCAGGGGCGCAGGCAAATCTAGCAGTGTATTTTGCAATGCTTGAACCGGGAGATAAGATTCTCGGTATGAACCTT
+
AAAA/E/EEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEAEA<EEE//EEAEE<EEAE/EEE/AAE/AEEEE6AEEE<AAEEEEE<AA6E/E//EE/A/E/6</EEAEEEA<<<AA//AEEE/AEE<A<

```

- **Paired-End Data Processing**

```
#Execution Command Examples

from MetaSAG import BarcodeDeal as bcd

#test_R1.fastq and test_R2.fastq
inputFastq = ["/DATA/test_R1.fastq", "/DATA/test_R2.fastq"] # Paired-end fastq files for analysis

CellBarn = Target_Path + "Barn/CellBarn_pair/" # Destination subdirectory for processed results (Ensure trailing "/" is present).

bcd.SAGSplit(inputFastq,CellBarn)


```

```
#Result

{'W1': 774, 'Len': 0, 'BC': 1533}
SAGSplit took 2.3399 seconds to execute.

```





## Func 2：trim(inputFastqDir,trimDir)
- **Function Description:**
Execute the **trimmomatic.jar**,The trimmomatic.jar tool is invoked to perform adapter trimming on sequencing reads for each cell in the inputFastqDir directory.

- **Required Parameters:**
```
inputFastqDir      Directory path containing droplet sequencing files
                   • Single-end FASTQ: Filename suffix must be .fastq (e.g., sample.fastq)
                   • Paired-end FASTQ: Filename suffixes must be _R1.fastq and _R2.fastq (e.g., sample_R1.fastq, sample_R2.fastq)

trimDir            Output directory for processed files after adapter trimming

```

- **Optional Parameters:**
```
ReadsEnd        Specifies whether input droplet sequencing files are single-end or paired-end
                Default: Single-end (ReadsEnd='Single')，


ILLUMINACLIP    --  Default parameter (JAR config): 'TruSeq3-PE.fa:2:30:10:3:TRUE'   

LEADING         --  Default parameter (JAR config): 25

TRAILING        --  Default parameter (JAR config): 3

SLIDINGWINDOW   --  Default parameter (JAR config): '4:20'

MINLEN          --  Default parameter (JAR config): 30

threads         --  Default parameter (JAR config): 12
```


Single-End Data Processing
```
#Execution Command Examples

from MetaSAG import BarcodeDeal as bcd

inputCellBarn = Target_Path + “Barn/CellBarn_single/”

trimBarn = Target_Path + “Barn/CellTrim_single/”

bcd.trim(inputCellBarn,trimBarn)

```


```
#Result

trim took 420.1405 seconds to execute.

```


Paired-End Data Processing
```
#Execution Command Examples

from MetaSAG import BarcodeDeal as bcd

inputCellBarn = Target_Path + “Barn/CellBarn_pair/”

trimBarn = Target_Path + “Barn/CellTrim_pair/”

bcd.trim(inputCellBarn,trimBarn,ReadsEnd='Pair')

```


```
#Result

trim took 663.4608 seconds to execute.

```