import os
import sys
import time
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF
import warnings

'''
脚本注释：
输入：
1，toy.fastq文件的目录位置 /xxx/xxx/xxx/toy.fastq
    toy.fastq:
    @Sam0516_10000:A00583:1296:HVHTKDSX7:2:1101:7871:16814 1:N:0:CGGAGCCT|GTGCATATTCGAGCACCGAATC
    CCTGTTACTTTACTAGGAACTACTCGTCACTTCTCTGGTGAAGAAGCTTCGGATTATCTGGTTGCAAAGAGTGATTGCC
    +
    FFFFFFFFFFFFFFF,FFFFFFFFFFFFF:FFFFFFFFFFFFFFF::FFFFFFFFFF:FF,FFF:FFFFFFFFFFFF,F
    @Sam0516_10000:A00583:1296:HVHTKDSX7:2:1101:25789:31657 1:N:0:CGGAGCCT|GTGCATATTCGAGCACCGAATC
    CATCATTCATTGGCCTAACGGAAAGTTATCTGATGGCCGAGTCAACACCAGCACCTTCCAATTGGTTAACGTTGCTCCA
    +
    FFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFF
    ... ...

2，BCFilterResult结果文件夹生成的目录位置 /xxx/xxx/xxx/BCFilterResult/ 以及前缀名称

输出：
BCFilterResult结果文件夹生成的目录位置 /xxx/xxx/xxx/BCFilterResult/ 以及前缀名称eg. toy
1，toy_bcread:
    cell	num
    500000	332
    500001	131
    500002	141
    500003	206
 2，1stDerivate.pdf
 3，2ndDerivate.pdf
 4，CDF_BCReads.pdf
 5，CellStack.pdf
 6，toy_BCFilter
    cell	num
    500019	88281
    500027	584
    500043	96147
    500048	14695


'''



'''
脚本注释
整合原BCFilter.py的内容。用户仅输入fastq文件的位置 + 结果文件位置，
返回值为min_reads，以及保留的细胞列表


'''



class BCFilter():
    def __init__(self,inputfastq,outputdir):
        self.inputfastq = inputfastq
        self.outputdir = outputdir

        # 判断inputfastq是单端还是双端
        if type(inputfastq) == str and inputfastq.endswith('.fastq'):
            self.ReadsEnd = 'Single'
            self.inputfastq = inputfastq
        elif type(inputfastq) == list:
            if len(inputfastq) == 2 and inputfastq[0].endswith('_R1.fastq') and inputfastq[1].endswith('_R2.fastq'):
                self.ReadsEnd = 'Pair'
                self.inputfastq1 = inputfastq[0]
                self.inputfastq2 = inputfastq[1]

            elif len(inputfastq) == 2 and inputfastq[0].endswith('_R2.fastq') and inputfastq[1].endswith('_R1.fastq'):
                self.ReadsEnd = 'Pair'
                self.inputfastq1 = inputfastq[1]
                self.inputfastq2 = inputfastq[0]

            elif len(inputfastq) == 1 and inputfastq[0].endswith('.fastq'):
                self.ReadsEnd = 'Single'
                self.inputfastq = inputfastq
            else:
                warnings.warn("传入的inputfastq不符合格式要求", UserWarning)


        else:
            warnings.warn("传入的inputfastq不符合格式要求",UserWarning)


    def timeit(func):
        def wrapper(*args,**kwargs):
            start_time = time.time()
            result = func(*args,**kwargs)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
            return result
        return wrapper




    @timeit
    def CellCountStatistic(self):
        ######Step1,统计prefix_bcread######

        #创建结果文件夹

        outputDir = self.outputdir
        ReadsEnd = self.ReadsEnd

        if ReadsEnd == 'Single':
            inputFastq = self.inputfastq
        else:
            inputFastq = self.inputfastq1


        if not os.path.exists(outputDir):
            os.makedirs(outputDir,exist_ok=True)


        BC_Count={}
        #读取inputFastq
        with open(inputFastq,'r') as f:
            while True:
                line = f.readline()
                if len(line) == 0:
                    break

                if line.startswith('@'):
                    line = line.replace('@','')
                    cell = line.split(':')[0]
                    if cell not in BC_Count.keys():
                        BC_Count[cell]=1
                    else:
                        BC_Count[cell]+=1

        BC_Count=pd.DataFrame(list(BC_Count.items()))
        BC_Count.columns=['cell','num']
        BC_Count.to_csv(os.path.join(outputDir,'bcread.txt'),sep='\t',index=False)
        self.BC_Count=BC_Count


    def getMinReads(self):
        outputDir = self.outputdir
        ######Step2,累积分布求导绘图######
        #BC_Count = pd.read_csv(BC_Count,sep='\t',header=0)
        BC_Count = self.BC_Count
        BC_Count_Sorted = BC_Count.sort_values('num')
        total_sum=BC_Count_Sorted['num'].sum()
        BC_Count_Sorted['consum'] = np.cumsum(BC_Count_Sorted['num'])
        BC_Count_Sorted['conprop'] = BC_Count_Sorted['consum']/total_sum

        x=[i for i in BC_Count_Sorted['num'] if i !=0]
        ecdf = ECDF(x)
        cell_readsnum=ecdf.x
        cumprop=ecdf.y

        plt.figure(figsize=(8, 8))
        plt.step(cell_readsnum, cumprop, where='post', label='ECDF')
        plt.xlabel('Cell_ReadsNum')
        plt.ylabel('CumProp')
        plt.title('Empirical Cumulative Distribution Function (ECDF)')
        plt.legend(loc='best')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'CDF_BCReads.pdf'))
        plt.close()
        #plt.show()

        ######Step3,最大优化算法求min_reads######
        reads_prop=pd.DataFrame(columns=['reads_num', 'prop'])

        for num in BC_Count_Sorted['num'].unique():
            temp1=num
            temp2=float(ecdf(num))
            row=[temp1,temp2]
            reads_prop=reads_prop._append(pd.Series(row,index=reads_prop.columns),ignore_index=True)


        reads_prop['n']=reads_prop['prop']-(reads_prop['reads_num'])/(len(reads_prop))
        flag = reads_prop.loc[reads_prop['n'] != np.inf, 'n'].max()
        min_reads=reads_prop.loc[(reads_prop['n'] == flag) & (reads_prop['n'] != np.inf),'reads_num'].values[0]



        ######Step4,统计一阶导,二阶导并绘图######
        #统计一阶导
        d1=[]
        for i in range(len(reads_prop)):
            if i==0:
                d=reads_prop.iloc[0]['prop']/1
            else:
                d=(reads_prop.iloc[i]['prop']-reads_prop.iloc[i-1]['prop'])/(reads_prop.iloc[i]['reads_num']-reads_prop.iloc[i-1]['reads_num'])

            d1+=[d]

        #d1PlotData=pd.DataFrame({'reads_num':reads_prop['reads_num'],'d1':d1})
        x=reads_prop['reads_num']
        y=d1
        plt.scatter(x,y,s=1)
        plt.xlim(-1000,5000)
        plt.xlabel('reads_num')
        plt.ylabel('D1')
        plt.title('Scatter Plot of reads number and D1')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'1stDerivate.pdf'))
        plt.close()
        #plt.show()
        #统计二阶导

        d2=[]
        for i in range(len(d1)):
            if i==0:
                d=d1[0]
            else:
                d=(d1[i]-d1[i-1])/1
            d2+=[d]


        x=reads_prop['reads_num']
        y=d2
        plt.scatter(x,y,s=1)
        plt.xlim(-1000,5000)
        plt.xlabel('reads_num')
        plt.ylabel('D1')
        plt.title('Scatter Plot of reads number and D2')
        plt.grid(True)
        plt.savefig(os.path.join(outputDir,'2ndDerivate.pdf'))
        plt.close()
        #plt.show()

        self.min_reads = min_reads

    def BCMinReads(self):
        BC_Count = self.BC_Count
        min_reads = self.min_reads
        outputDir = self.outputdir

        ######Step5,统计过滤得到的细胞编号######
        #BC_Count = pd.read_csv(BC_Count, sep='\t', header=0)
        BC_Count_Filter = BC_Count[BC_Count['num']>=min_reads]
        BC_Count_Filter.to_csv(os.path.join(outputDir,'BC_Count_Filter.txt'),sep='\t',index=False)
        self.BC_Count_Filter = BC_Count_Filter

    '''
        if inputFastq is not None:
            if type(inputFastq) == str:
                inputFastq = [inputFastq]
    
            for ifq in inputFastq:
                if ifq.endswith('_R1.fastq'):
                    BCFilterFile=os.path.join(outputDir,'BCFilter_R1.fastq')
                elif ifq.endswith('_R2.fastq'):
                    BCFilterFile = os.path.join(outputDir, 'BCFilter_R2.fastq')
                elif ifq.endswith('.fastq'):
                    BCFilterFile = os.path.join(outputDir, 'BCFilter.fastq')
    
                with open(ifq,'r') as f:
                    while True:
                        line = f.readline()
                        if len(line) == 0:
                            break
                        line = line.strip('\n')
    
                        if line.startswith('@'):
                            line1 = line.replace('@', '')
                            cell = line1.split(':')[0]
    
                        if cell in BC_Count_Filter['cell'].values:
                            flag = True
                        else:
                            flag = False
    
                        if flag:
                            with open(BCFilterFile,'a') as f2:
                                f2.write(line+'\n')
        
    '''



    def Sam2Cell(self):
        inputfastq=self.inputfastq
        CellBarnDir=os.path.join(self.outputdir,'CellBarnDir')
        #创建CellBarnDir的结果目录
        if not os.path.exists(CellBarnDir):
            os.makedirs(CellBarnDir,exist_ok=True)


        for ifq in self.inputfastq:
            if ifq.endswith('_R1.fastq'):
                suffix = '_R1.fastq'
            elif ifq.endswith('_R2.fastq'):
                suffix = '_R2.fastq'
            elif ifq.endswith('.fastq'):
                suffix = '.fastq'

            with open(ifq,'r') as f:
                while True:
                    line = f.readline()
                    if len(line) == 0:
                        break
                    line = line.strip('\n')

                    if line.startswith('@'):
                        line1 = line.replace('@','')
                        cell = line1.split(':')[0]

                    with open(os.path.join(CellBarnDir,cell)+suffix,'a') as output:
                        output.write(line+'\n')




def SamsBCFilterStack(Sam,SavedCell,AllCell,outputDir):
    #AllCell = [3132, 4557, 2321, 3433, 4500]
    #SavedCell = [2000,3000,1734,2222,2000]
    #Sam=['S5','S6','S7','S8','S10']
    fig, ax = plt.subplots()
    ax.bar(Sam,AllCell,label='Filtered Cells')
    ax.bar(Sam, SavedCell, label='Saved Cells')
    ax.set_title('Stacked Bar Chart')
    ax.set_xlabel('Sam')
    ax.set_ylabel('Cell Number')
    ax.legend()

    # 显示图形
    plt.savefig(os.path.join(outputDir,'CellStack.pdf'))
    plt.close()
    #plt.show()



'''
if __name__ == '__main__':
    from zszshhh import BCFilterClass
    import pandas as pd
    x=BCFilterClass.BCFilter(['D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\BCFilter\\input\\toy_R2.fastq','D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\BCFilter\\input\\toy_R1.fastq'],'D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\BCFilter\\result')    BC_Count = pd.read_csv('D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\BCFilter\\input\\S5_bcread',sep='\t',header=0)
    x.ReadsEnd
    x.inputfastq1
    x.BC_Count=BC_Count
    x.getMinReads()
    x.BCMinReads()
    x.Sam2Cell()
'''
