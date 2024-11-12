# Step 1. Distribute the reads in the sample to a file of individual droplets.

## Func 1：SAGSplit(inputFastq,CellBarn)

- **函数功能：**

SAGSplit()输入的样本测序文件（fastq格式），根据barcode序列信息，将每条reads分别写入到相应的液滴文件中。

reads的barcode识别由内置的FindBarcode函数（或用户自己提供）完成。


- **必选参数：**
```
inputFastq  --  样本的短reads测序文件位置
                如果是单端fastq文件，文件名必须以.fastq结尾；
                如果是双端fastq文件，文件名必须以_R1.fastq或_R2.fastq结尾，并以列表格式给出。

CellBarn    --  分割样本后存储液滴文件的结果路径
```

- **可选参数：**
```
FindBarcode --  要求输入的是一个用户自己写的函数，该函数只输入一个reads字符串，要求返回一个字符串，该字符串在正常匹配时返回reads所属的Barcode字符串，或者未匹配的原因字符串eg,"Erro","Len","W1","BC"总之能返回一个字符串。<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                默认为脚本中定义的函数FindBarcode()

warning --  FindBarcode函数中设定某些reads无法正常匹配Barcode时自编函数返回字符串的情况
            默认为['Len','W1','BC']

filterWarning   --  如果reads对应的Barcode字符串在warning中，这样的reads写入到单液滴文件中
                    默认为True
```


- **对于单端数据**
```
#输入文件示例
test.fastq:

@NB501288_516_HC2NTBGXB:1:11101:18376:1042#CGAGGCTG/1
AGAGGTNGGAGTGATTGCTTGTGACGCCTTTGCCTCACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGTCCTTAACCATCCTTGAATACCTCGCTTGCTATTTTTTGTGCTTCTTTCCTCAGATATTGTGCCGTCTCATACATGAATGGTCTTT
+
AAAAAE#EEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAEEEEEAEEEEEEEEEE6AEEEEAEEEAEEEEEE/EE/EEEEAEAAEE/<EEE<AE</EEE<AAEE6EEAAEEE<EAEEEEEAEAEAE/A/EEAA</6/<EEE/EEEE</</E6AA/
@NB501288_516_HC2NTBGXB:1:11101:20037:1043#CGAGGCTG/1
GTTTGTNTGAGTGATTGCTTGTGACGCCTTGTTTGTTTTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGATAAATACGTATAGTACGATCAAAAACGCAAGAATATATCCGATCGCCCCAAGCGCCGGAATGCCGAGTATTTTCGGACTCATATC
+
AAAAAE#EEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE6/EEEEEEAEEAAEEE/EEEEEEEEEEEEEEEEEEEEEEEEA<AEEEEAAEEE<A/EE<<<<<AE/E/A<A<A<EEA66AAEEEEEAE<<<E<A<A6AA/

··· ···
```


```
#执行代码
import BarcodeDeal as bcd
inputFastq = 'zszshhh/testData/BarcodeDeal/input/test.fastq'
CellBarn = './testData/BarcodeDeal/result/CellBarn_single'
bcd.SAGSplit(inputFastq,CellBarn)

```


```
#结果文件示例
#318.fastq
#该文件中每条reads的标头以@Barcode_ID:开头

@318:NB501288_516_HC2NTBGXB:1:11101:8220:5572#CGAGGCTG/1
GTTTGTTTGAGTGATTGCTTGTGACGCCTTCCTGACACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGCTTGTATACAATATGCTTATAGTATACTCATATTTTCCTTAAAAATCAATATTTTATCTCACGATTTTAAATCTGAATTTTCCATTT
+
AAAAA/EEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEEEEEEEEEEAEEEEEE6AEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE<AEAEAEEAE<<EEAEEEEE6AEEEE/A/E<<<///AEE66AEEEEAEAEEAE/A</<</A
@318:NB501288_516_HC2NTBGXB:1:11101:23509:16472#CGAGGCTG/1
GTTTGTTTGAGTGATTGCTTGTGACGCCTTCCTGACACTCGTCGGCAGCGTCAGATGTCTATAAGAGACAGGTACAGCCGCATTCAGGGGCGCAGGCAAATCTAGCAGTGTATTTTGCAATGCTTGAACCGGGAGATAAGATTCTCGGTATGAACCTT
+
AAAA/E/EEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEAEA<EEE//EEAEE<EEAE/EEE/AAE/AEEEE6AEEE<AAEEEEE<AA6E/E//EE/A/E/6</EEAEEEA<<<AA//AEEE/AEE<A<

```

- **对于双端数据**

```
#执行代码
import BarcodeDeal as bcd
inputFastq = ['zszshhh/testData/BarcodeDeal/input/test_R1.fastq','zszshhh/testData/BarcodeDeal/input/test_R2.fastq']
CellBarn = './testData/BarcodeDeal/result/CellBarn_pair'
bcd.SAGSplit(inputFastq,CellBarn)

```



## Func 2：trim(inputFastqDir,trimDir)
- **函数功能：**
调用jar包 **trimmomatic.jar**,对inputFastqDir目录下每个细胞中的reads进行去接头的修剪步骤。

- **必选参数：**
```
inputFastqDir   --  每个液滴测序文件存放的路径<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
如果是单端fastq文件，文件名必须以.fastq结尾；<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
如果是双端fastq文件，文件名必须以_R1.fastq或_R2.fastq结尾。

trimDir -- 接头去除后的结果文件存放位置
```

- **可选参数：**
```
ReadsEnd    --  输入液滴测序文件是单端还是双端<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
默认为单端，ReadsEnd='Single'；<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
如果是双端，修改ReadsEnd='Pair'.


ILLUMINACLIP    --  jar包参数,默认为'TruSeq3-PE.fa:2:30:10:3:TRUE'   

LEADING --  jar包参数,默认为25

TRAILING --  jar包参数,默认为3

SLIDINGWINDOW --  jar包参数,默认为'4:20'

MINLEN --  jar包参数,默认为30

threads --  jar包参数,默认为12
```