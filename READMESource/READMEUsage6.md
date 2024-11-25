# MetaSAG Usage 
## Step 6. Species to Strain resolved genomes.

## Class1：SingleBin(BinDir, ResultDir, SpeciesName)
- **类功能：**

对单个物种水平箱中的细胞分簇为菌株水平。

- **必选参数：**
```
BinDir      --      物种箱路径
                    该箱中要求包含每个细胞的fastq测序文件，以及箱组装基因组fasta文件。    

ResultDir   --      分为菌株簇后的结果路径

SpeciesName --      该箱组装基因组文件名(去除.fasta后缀)，同时为结果文件名前缀。

```

## Func 1：SingleBinPrepare()

- **函数功能：**

对输入箱中的细胞进行变异召回，根据SNP对细胞分簇绘图。

- **可选参数：**

```

bcftools        --      调用bcftools工具路径。
                        默认为None

snap_aligner    --      调用snap_aligner工具路径。
                        默认为None

env             --      bcftools/snap_aligner运行需要的conda环境
                        默认为None

ReadsEnd        --      输入液滴测序文件是单端还是双端
                        默认为单端，ReadsEnd='Single'
                        如果是双端，修改ReadsEnd='Pair'.

```




- **结果：**

![SNP_Pheatmap](SNP_Pheatmap.png)
![SNP_Tree](SNP_Tree.png)
![SNP_Umap_Prep](SNP_Umap_prep.png)







## Func 2：SingleBinSplit(ClusterNum)

- **函数功能：**

根据SingleBinPrepare函数得到的图片结果，判断该箱应该分为几簇。

- **必选参数：**

ClusterNum      --      该箱被分为合适的簇数

- **结果：**

![SNP_Umap](SNP_Umap.png)


Eg. StrainCell.txt

| Cluster |     Cell      |
|:-------:|:-------------:|
|    0    | Sam1025_10073 |
|    0    | Sam1025_1115  |
|    0    | Sam1025_11610 |
|    1    | Sam1025_11666 |
|    1    | Sam1025_12999 |
|    1    | Sam1025_4770  |
|   ...   |      ...      |

```
#执行代码示例

from MetaSAG import SNPStrain as snp


#分箱预处理

BinDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/input/SingleBin/SGB6796' #1.5Gb

ResultDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/result/SingleBin/SGB6796'

SpeciesName='SGB6796'  

SingleBin=snp.SingleBin(BinDir,ResultDir,SpeciesName)

snap_aligner = '/data_alluser/singleCellMicrobiome/dmy_test/tools/SNAP/snap-aligner'

bcftools = '/data_alluser/singleCellMicrobiome/dmy_test/tools/bcftools/bcftools-1.18/bcftools'

SingleBin.SingleBinPrepare(bcftools=bcftools,snap_aligner=snap_aligner,ReadsEnd='Pair')
#SingleBinPrepare took 216.9003 seconds to execute.


#查询去除SNP/细胞的个数

SingleBin.DropSNPNum  # 12057

SingleBin.AllSNPNum  #  22728

SingleBin.AllCellNum  #  69

SingleBin.DropCellNum  #  3




#分簇
#根据SingleBinPrepare()得到的umap图确定分割簇的个数

SingleBin.SingleBinSplit(2)
#SingleBinSplit took 1.2961 seconds to execute.
```


## Class2：AllBin(FastaDir, FastqDir, CellAnno, ResultDir)
- **类功能：**

对单个物种水平箱中的细胞分簇为菌株水平。

- **必选参数：**
```
FastaDir        --      所有箱组装基因组文件存放路径    

FastqDir        --      所有细胞fastq测序文件存放路径

CellAnno        --      每个细胞与箱名之间的对应关系

ResultDir       --      分箱结果路径

```

## Func 1：AllBinPrepare()

- **函数功能：**

对输入箱中的细胞进行变异召回，根据SNP对细胞分簇绘图。

- **可选参数：**
```

bcftools        --      调用bcftools工具路径。
                        默认为None

snap_aligner    --      调用snap_aligner工具路径。
                        默认为None

env             --      bcftools/snap_aligner运行需要的conda环境
                        默认为None

ReadsEnd        --      输入液滴测序文件是单端还是双端
                        默认为单端，ReadsEnd='Single'
                        如果是双端，修改ReadsEnd='Pair'.

```


Eg. CellAnno.txt

| Cluster  |     Cell      |
|:--------:|:-------------:|
| SGB10068 | Sam1025_10158 |
| SGB10068 | Sam1025_10251 |
| SGB10068 | Sam1025_10392 |
| SGB4577  | Sam1102_7184  |
| SGB4577  | Sam1102_7261  |
|   ...    |      ...      |


## Func 2：AllBinSplit()

- **函数功能：**

根据AllBinPrepare函数得到的图片结果，判断每个箱应该分为几簇。
将合适的分簇数填写到ResultDir/BinClusterAnno.txt文件中（或者在运行函数前修改AllBin对象的BinClusterAnno,即文件路径)


Eg. ResultDir/BinClusterAnno.txt

|         BinPrepareDir          | SpeciesName | ClusterNum |
|:------------------------------:|:-----------:|:----------:|
| ./result/AllBin/SGB10068Result |  SGB10068   |     4      |
| ./result/AllBin/SGB5111Result  |   SGB5111   |     2      |
|              ...               |     ...     |    ...     |


```
#执行代码示例

from MetaSAG import SNPStrain as snp


#分箱预处理

#fasta + fastq 14Gb
fastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/input/AllBin/fasta'

fastqDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/input/AllBin/fastq'

cellAnno='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/input/AllBin/SNPCluster.txt'

resultDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/SNPStrain/result/AllBin'

AllBin=snp.AllBin(fastaDir,fastqDir,cellAnno,resultDir)

snap_aligner = '/data_alluser/singleCellMicrobiome/dmy_test/tools/SNAP/snap-aligner'

bcftools = '/data_alluser/singleCellMicrobiome/dmy_test/tools/bcftools/bcftools-1.18/bcftools'

AllBin.AllBinPrepare(bcftools=bcftools,snap_aligner=snap_aligner,ReadsEnd='Pair')
#AllBinPrepare took 3736.5503 seconds to execute.

#分箱
#注意分箱前根据结果图片,在ResultDir/BinClusterAnno.txt文件中填写每个箱合适的分簇数

AllBin.AllBinSplit()
#AllBinSplit took 43.4115 seconds to execute


```