# MetaSAG Usage 
## Step 9. Droplet clustering of potentially unknown species

## Func 1：ClusterSAG(inputFastq,outputDir)
- **函数功能：**

对筛选出的潜在未知物种进行聚类（第一次聚类）。

- **必选参数：**
```

inputFastq      --      潜在未知物种液滴fastq文件存放路径

outputDir       --      聚类结果文件存放路径

```

- **可选参数：**
```
ReadsEnd        --      输入液滴测序文件是单端还是双端
                        默认为单端，ReadsEnd='Single'
                        如果是双端，修改ReadsEnd='Pair'.
                        
SpadesEnv       --      Spades.py运行依赖的conda环境
                        默认为None
                        
SourmashEnv     --      Sourmash运行依赖的conda环境
                        默认为None                   

```



## Func 2：ClusterBin(OldRoundFold,NewRoundFold)

- **函数功能：**

根据用户需要，在之前聚类的结果上进行再一次聚类。


- **必选参数：**
```

OldRoundFold    --      上一次聚类结果存放的根路径（函数输入数据）

outputDir       --      再次聚类的结果文件存放路径

```

- **可选参数：**
```
ReadsEnd        --      输入液滴测序文件是单端还是双端
                        默认为单端，ReadsEnd='Single'
                        如果是双端，修改ReadsEnd='Pair'.
                        
SpadesEnv       --      Spades.py运行依赖的conda环境
                        默认为None
                        
SourmashEnv     --      Sourmash运行依赖的conda环境
                        默认为None                   

```



```

#执行代码示例

from MetaSAG import UnknownSAG as usag

fastqDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/MetaPhlAnAsign/result/UnknownCell/fastq' #846M

resultDir='/data_alluser/singleCellMicrobiome/dmy_test/gj/MetaPhIAn4_1/PyPack/PyPackData2/testData/MetaPhlAnAsign/result/UnknownCell/testResult'

usag.ClusterSAG(fastqDir,resultDir,SourmashEnv='sourmash') #对未知物种液滴进行初次聚类

Round1Dir=os.path.join(resultDir,'Round1')

Round2Dir=os.path.join(resultDir,'Round2')

usag.ClusterBin(Round1Dir,Round2Dir) #根据初次聚类结果对未知物种液滴的箱再次聚类

```