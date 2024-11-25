# MetaSAG Usage 
## Step 5. Build phylogenetic tree.


## Func 1：BuildTree(FastaDir,TreeTemp,env=None)

- **函数功能：**

调用ANVI'O,对输入目录下的基因组文件构建进化树.


- **必选参数：**
```
FastaDir        --      输入基因组文件路径
                        注意必须是绝对路径。

TreeTemp        --      树文件结果路径



```


- **可选参数：**
```
env         --      ANVI'O运行的conda环境。
                    默认为None。
```

- **结果：**

phylogenomic-tree_Bacteria_71_ribosomal6.txt




## Func 2：itolPlot(BinAnno,Anno)

- **函数功能：**

将提供的基因组信息整理为适合itol网页绘制进化树的数据格式。


- **必选参数：**
```

BinAnno     --      基因组文件信息

Anno        --      itol绘图格式结果文件存放路径


```

Eg. BinAnno (不提供门水平对应颜色)

|   Bin   | CellNum |     Phylum     |
|:-------:|:-------:|:--------------:|
| genome1 |   91    |   Firmicutes   |
| genome2 |   44    | Actinobacteria |
| genome3 |   18    | Actinobacteria |
|   ...   |   ...   |      ...       |



Eg. BinAnno (提供门水平对应颜色)

|   Bin   | CellNum |     Phylum     |  Color  |
|:-------:|:-------:|:--------------:|:-------:|
| genome1 |   91    |   Firmicutes   | #705e78 |
| genome2 |   44    | Actinobacteria | #fea443 |
| genome3 |   18    | Actinobacteria | #fea443 |
|   ...   |   ...   |      ...       |   ...   |

- **结果：**

![Tree](Tree.png)



```
#执行代码示例

from MetaSAG import Tree as tree

#构建树

FastaDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/Tree/input' #292Mb

TreeTemp='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/Tree/result'

tree.BuildTree(FastaDir,TreeTemp,env='anvio-7.1')
#BuildTree took 8312.0662 seconds to execute.



#整理itol网页树绘图文件

BinAnno='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/Tree/input/BinAnno'

AnnoResult='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/Tree/result/Anno'

tree.itolPlot(BinAnno,AnnoResult)


```