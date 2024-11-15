# MetaSAG Usage 
## Step 4. Quality control and integration annotation of assembled genome.

## Func 1：FastaQC1(FastaDir,FastaOut)

- **函数功能：**
去除Fasta路径下所有fasta文件中过短的contig.

- **必选参数：**
```
FastaDir        --      未校正fasta文件存放路径
                        要求文件以.fasta结尾

FastaOut        --      校正后fasta文件的结果路径

```````
- **可选参数：**
```
minlen          --      fasta文件中contig长度的最小阈值，
                        默认为500bp。

```


## Func 2：FastaQC2(FastaDir,FastaOut,CheckmFile)

- **函数功能：**
结合checkm结果，将校正后基因组fasta文件归类为合格/需要再校正/放弃。

- **必选参数：**
```
FastaDir        --      校正后fasta文件存放路径

FastaOut        --      fasta文件分类结果路径

CheckmFile      --      校正后fasta文件Checkm结果表格

```````
- **可选参数：**
```
FastgDir        --      如果能够提供fastg文件夹，将相应fastg文件复制到Bandage文件夹中； 如果没有提供fastg文件夹，则只将相应fasta文件放到Bandage文件夹中
                        默认为None。

```


## Func 3：Summary(FastaSGB,CheckmFile,GTDBFile)

- **函数功能：**

将用户提供的组装基因组的MetaPhlAnN4分类文件(FastaSGB),Checkm质量检查文件(CheckmFile),GTDB-TK注释文件(GTDBFile)的结果整合到一起。

- **必选参数：**
```
FastaSGB        --      组装基因组经MetaPhlAnN4注释的分类文件路径

CheckmFile      --      组装基因组的Checkm质量检查文件

GTDBFile        --      组装基因组的GTDB-TK注释文件

```


- **可选参数：**
```
outputSummary   --      注释结果整合结果文件路径
                        默认为None,注释结果形式为DataFrame对象。
```



Eg. FastaSGB(第一列为箱基因组文件名，要求去除后缀.fasta;第二列为MetaPhlan4分类SGB编号)

|   Bin   |   SGB    |
|:-------:|:--------:|
| genome1 | SGB1024  |
| genome2 | SGB28348 |
| genome3 | SGB4348  |
|   ...   |   ...    |


- **结果：**

Eg. Summary.txt

|   Bin   |   SGB    |                                                         SGBName                                                          | Completeness | Contamination | closest_reference |                                                     closest_taxonomy                                                     |closest_ani|ContamFlag|CompleteFlag|QC|
|:-------:|:--------:|:------------------------------------------------------------------------------------------------------------------------:|:------------:|:-------------:|:-----------------:|:------------------------------------------------------------------------------------------------------------------------:|:---:|:---:|:---:|:---:|
| genome1 | SGB1024  |              k__Bacteria\|p__Bacteroidetes\|c__CFGB343\|o__OFGB343\|f__FGB343\|g__GGB781\|s__GGB781_SGB1024              |    96.02     |     13.66     |  GCA_000434395.1  |                d__Bacteria;p__Firmicutes;c__Bacilli;o__RFN20;f__CAG-288;g__CAG-568;s__CAG-568 sp000434395                |97.29|fail|pass|Low|
| genome2 | SGB28348 | k__Bacteria\|p__Bacteroidetes\|c__Cytophagia\|o__Cytophagales\|f__Cytophagaceae\|g__Aquirufa\|s__Aquirufa_antheringensis |    86.21     |     71.32     |        no         |                                                            no                                                            |no|fail|pass|Low|
| genome3 | SGB4348  |            k__Bacteria\|p__Firmicutes\|c__CFGB4806\|o__OFGB4806\|f__FGB4806\|g__GGB51647\|s__GGB51647_SGB4348            |    98.66     |     7.31      |  GCA_000432435.1  | d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Fimenecus;s__Fimenecus sp000432435 |98.92|pass|pass|Medium|
|   ...   |   ...    |                                                           ...                                                            |     ...      |      ...      |        ...        |                                                           ...                                                            |...|...|...|...|




```
#执行代码示例

from MetaSAG import BinQCAnno as bqa



#去除长度<500bp的contig

FastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC1Fasta'

FastaOut='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/result/QC1Fasta'

bqa.FastaQC1(FastaDir,FastaOut)



#根据checkm文件判断基因组文件的组装质量。判断哪些组装基因组可以进一步校正

FastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/fasta'

FastgDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/fastg'

CheckmFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/QC2Fasta/qa_result'

FastaOut='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/result/QC2Fasta'

bqa.FastaQC2(FastaDir,FastaOut,CheckmFile,FastgDir=FastgDir) 



#整合MetaPhlAnN4注释分类,Checkm质量评估以及GTDBTK注释分类结果。

FastaSGB='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/FastaSGB'

CheckmFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/qa_result'

GTDBFile='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/gtdbtk.bac120.summary.tsv'

outputSummary='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/BinQCAnno/input/summary/summary.txt'

summary=bqa.Summary(FastaSGB,CheckmFile,GTDBFile)




```

