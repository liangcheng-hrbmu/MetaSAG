# MetaSAG Usage 
## Step 7. Horizontal Gene Transfer.

## Class：CellHGT(FastaDir,outputDir)
- **类功能：**

找出两两基因组文件(物种水平)间的水平基因转移。

- **必选参数：**
```
FastaDir        --      输入基因组文件路径

outputDir       --      结果路径

```

## Func 1：SpeciesHGT()

- **函数功能：**
计算两两基因组文件(物种水平)间的水平基因转移，返回HGT.fasta文件
要求FastaDir下的基因组文件命名为物种名.fasta(或物种编号.fasta),不允许存在下划线。


Eg. HGT.fasta

```

>PAIR_物种名1_ContigIDTO物种名2_ContigID
CCACCATGTATGACTGGCTTGCCACGATT...
>PAIR_SGB4910_36TOSGB15286_37
GTCTATTGATGAGCAAGGACTGAGCAGTG...
...

```


## Func 2：HGTSpeciesPlot(TreeAnno)

- **函数功能：**

将提供的树节点信息整理为适合itol网页绘制进化树的数据格式。

Eg. TreeAnno

|   File    | SpeciesID |     Phylum     | Phylum_Color |
|:---------:|:---------:|:--------------:|:------------:|
| SGB15452  | SGB15452  | Proteobacteria |   #fea443    |
|  SGB9262  |  SGB9262  | Proteobacteria |   #fea443    |
|  SGB6019  |  SGB6019  |  Fusobacteria  |   #b46b6b    |
| SGB102029 | SGB102029 |   Firmicutes   |   #705E78    |
| SGB14991  | SGB14991  |   Firmicutes   |   #705E78    |
|    ...    |    ...    |      ...       |     ...      |


![HGT_Tree](HGT_Tree.png)




## Func 3：HGTSpeciesAnno(TreeAnno)

- **函数功能：**

对HGT.fasta文件注释基因并进行相似性聚类,统计每簇基因涉及到的物种。

- **必选参数：**
```
TreeAnno        --      树节点注释文件

```

- **可选参数：**

```
prokka_env      --      prokka运行的conda环境
                        默认为None

cdhit_env       --      cdhit运行的conda环境
                        默认为None

emapper_env     --      Eggnog-mapper运行的conda环境
                        默认为None

emapper_DB      --      Eggnog-mapper参考数据库路径
                        默认为None

```


- **结果：**

Eg. ./SpeciesHGTResult/HGTSpeciesAnno/AnnoCDHit/CDHitCluster.txt

|  Cluster  |      Gene      |           Contig            |  Bin1   |  Bin2   |
|:---------:|:--------------:|:---------------------------:|:-------:|:-------:|
| Cluster 0 | HDNNIBCP_01849 |  PAIR_SGB4573_3TOSGB4826_1  | SGB4573 | SGB4826 |
| Cluster 0 | HDNNIBCP_02231 | PAIR_SGB4573_3TOSGB4874_121 | SGB4573 | SGB4874 |
| Cluster 1 | HDNNIBCP_04222 | PAIR_SGB4874_121TOSGB4826_1 | SGB4874 | SGB4826 |
|    ...    |      ...       |             ...             |   ...   |   ...   |



Eg. ./SpeciesHGTResult/HGTSpeciesAnno/AnnoCDHit/ClusterBin.txt

|  Cluster  |         Binlist         |
|:---------:|:-----------------------:|
| Cluster 0 | SGB4874,SGB4573,SGB4826 |
| Cluster 1 |     SGB4874,SGB4826     |
|    ...    |           ...           |




## Func 4：StrainHGT()

- **函数功能：**

计算两两基因组文件(菌株水平)间的水平基因转移，返回HGT.fasta文件
要求FastaDir下的基因组文件命名为物种名@菌株编号.fasta(或物种编号@菌株编号.fasta),不允许存在下划线。

Eg. HGT.fasta

```

>PAIR_物种名1@菌株号_ContigIDTO物种名2@菌株号_ContigID
CCACCATGTATGACTGGCTTGCCACGATT...
>PAIR_yw14@strain0_24TOyw90@strain0_14
GGTTCTTGTAGTTGTGGGCCTCGTCCACAAACAGCCGGT
...

```


## Func 5：HGTStrainAnno(TreeAnno)

- **函数功能：**

对HGT.fasta文件注释基因并进行相似性聚类,统计每簇基因涉及到的菌株。

- **必选参数：**
```
TreeAnno        --      树节点注释文件

```

- **可选参数：**

```
prokka_env      --      prokka运行的conda环境
                        默认为None

cdhit_env       --      cdhit运行的conda环境
                        默认为None

emapper_env     --      Eggnog-mapper运行的conda环境
                        默认为None

emapper_DB      --      Eggnog-mapper参考数据库路径
                        默认为None

```


```
#执行代码示例

from MetaSAG import CellHGT as hgt



#得到两两基因组(物种水平)之间的水平基因转移序列 HGT.fasta

fastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/CellHGT/input/testSpeciesFasta'

HGTTemp='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/CellHGT/result/SpeciesHGTResult'

obj=hgt.CellHGT(fastaDir,HGTTemp)

obj.SpeciesHGT()




##提供itol网页需要的绘图注释文件

TreeAnno='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/CellHGT/input/SpeciesAnno'

obj.HGTSpeciesPlot(TreeAnno)


#对水平基因转移Contig注释基因并根据相似性聚类统计。

obj.HGTSpeciesAnno(TreeAnno,prokka_env='prokka',cdhit_env='base',emapper_env='eggnog-mapper2',emapper_DB='/data_alluser/public/database/eggnogDB/')


```


