import math
import os
import pandas as pd
import numpy as np
import shutil
import subprocess

'''
脚本注释
本脚本要求输入一个fasta文件夹。
1. 提供每个fasta文件对应的SGB编号
2. 使用本Python包中自带的MetaPhlan4的SGB注释文件对每个fasta文件进行物种分类命名的整理
3. Checkm部分选择是否进行Bandage二次校正。如果进行二次校正，需要提供自己本来污染超标但是Bandage校正后质量达标的fasta文件夹


'''

PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
def removeShortContig(inputfastafile,outputfastafile,minlen=500):
    fasta={}
    with open(inputfastafile) as input1:
        while True:
            line = input1.readline()
            if len(line) == 0:
                break

            line = line.replace('\n', '')

            if line.startswith('>'):
                flag = line
                fasta[flag] = []
            else:
                fasta[flag] += [line]
    input1.close()

    # 长度大于等于500的contigs
    lt_500 = []
    for key in fasta.keys():
        l = int(key.split('_')[3])
        if l >= minlen:
            lt_500 += [key]

    # 覆盖多于两个标准差
    cov = []
    for node in lt_500:
        cov += [float(node.split('_')[5])]

    cov_log = [math.log(num) for num in cov]
    pj = np.mean(cov_log)
    bz = np.std(cov_log)
    threshold = pj - 2 * bz

    cov_save = []
    for i in range(len(lt_500)):
        if cov_log[i] > threshold:
            cov_save += [lt_500[i]]

    with open(outputfastafile, 'a') as output:
        for key, value in fasta.items():
            if key in cov_save:
                output.write(key + '\n')
                for line in value:
                    output.write(line + '\n')

    output.close()





def Summary(FastaSGB,CheckmFile,GTDBFile,outputSummary=None):

    mpa = pd.read_csv(os.path.join(PYTHONDIR,'mpa_vOct22_CHOCOPhlAnSGB_202403_species.txt'),sep='\t',header=None)
    mpa.columns=['SGB','SGBSpecies']

    FastaSGB = pd.read_csv(FastaSGB,sep='\t',header=0)
    FastaSGB.columns = ['Bin','SGB']
    Checkm = pd.read_csv(CheckmFile,header=0,sep='\t')
    GTDB = pd.read_csv(GTDBFile,header=0,sep='\t')

    GTDB = GTDB[['user_genome','closest_placement_reference','closest_placement_taxonomy','closest_placement_ani']]
    GTDB.columns = ['Bin','closest_reference','closest_taxonomy','closest_ani']
    Checkm = Checkm[['Bin Id','Completeness','Contamination']]
    Checkm.columns = ['Bin','Completeness','Contamination']

    SummaryResult = pd.merge(FastaSGB,mpa,on='SGB',how='left')
    SummaryResult = pd.merge(SummaryResult,Checkm,on='Bin',how='left')
    SummaryResult = pd.merge(SummaryResult,GTDB,on='Bin',how='left')

    #将没有匹配的所有单元格填充为'no'
    SummaryResult.fillna('no',inplace=True)

    SummaryResult['ContamFlag'] = None
    SummaryResult['CompleteFlag'] = None
    SummaryResult['QC'] = None
    for index, row in SummaryResult.iterrows():
        if row['Contamination'] != 'no':
            if row['Contamination'] > 10:
                SummaryResult.loc[index,'ContamFlag'] = 'fail'
            else:
                SummaryResult.loc[index, 'ContamFlag'] = 'pass'
        else:
            SummaryResult.loc[index, 'ContamFlag'] = 'no'


        if row['Completeness'] != 'no':
            if row['Completeness'] < 50:
                SummaryResult.loc[index,'CompleteFlag'] = 'fail'
            else:
                SummaryResult.loc[index, 'CompleteFlag'] = 'pass'
        else:
            SummaryResult.loc[index, 'CompleteFlag'] = 'no'


        if row['Contamination'] != 'no' and row['Completeness'] != 'no':
            if (row['Contamination'] <= 5) & (row['Completeness'] >= 90):
                SummaryResult.loc[index,'QC'] = 'High'
            elif (row['Contamination'] <= 10) & (row['Completeness'] >= 50):
                SummaryResult.loc[index, 'QC'] = 'Medium'
            else:
                SummaryResult.loc[index, 'QC'] = 'Low'
        else:
            SummaryResult.loc[index, 'QC'] = 'no'



    '''
    SummaryResult['ContamFlag'] = np.where(SummaryResult['Contamination'] > 10,'fail','pass')
    SummaryResult['CompleteFlag'] = np.where(SummaryResult['Completeness'] < 50,'fail','pass')

    conditions = [
        (SummaryResult['Contamination'] <= 5) & (SummaryResult['Completeness'] >= 90),
        (SummaryResult['Contamination'] <= 10) & (SummaryResult['Completeness'] >= 50)
    ]
    choices = ['High','Medium']
    SummaryResult['QC'] = np.select(conditions,choices,default='Low')
    '''

    if outputSummary is not None:
        SummaryResult.to_csv(outputSummary,sep='\t',index=False)

    return SummaryResult



def FastaQC1(FastaDir,FastaOut):
    #该函数只进行去除短Contig的操作
    #创建目录
    if not os.path.exists(FastaOut):
        os.makedirs(FastaOut,exist_ok=True)

    files = os.listdir(FastaDir)
    for filename in files:
        if filename.endswith('.fasta'):
            removeShortContig(os.path.join(FastaDir, filename),os.path.join(FastaOut, filename))



def FastaQC2(FastaDir,FastaOut,CheckmFile,FastgDir=None):
    #该函数
    #(1)对去除短contig的fasta文件根据Checkm结果文件，返回哪些文件是合格基因组[Pass]，哪些文件虽然不合格但是可以Bandage校正[Bandage]，哪些基因组文件由于完整度过小可以直接放弃[Abandon]
    #(2)如果能够提供fastg文件夹，将相应fastg文件复制到Bandage文件夹中； 如果没有提供fastg文件夹，则只将相应fasta文件放到Bandage文件夹中

    #创建目录
    Pass = os.path.join(FastaOut,'Pass')
    Bandage = os.path.join(FastaOut,'Bandage')
    Abandon = os.path.join(FastaOut,'Abandon')

    DirCreate = [FastaOut,Pass,Bandage,Abandon]
    for dir in DirCreate:
        if not os.path.exists(dir):
            os.makedirs(dir,exist_ok=True)

    Checkm = pd.read_csv(CheckmFile, header=0, sep='\t')
    Checkm = Checkm[['Bin Id', 'Completeness', 'Contamination']]
    Checkm.columns = ['Bin', 'Completeness', 'Contamination']

    for index,row in Checkm.iterrows():
        Bin = row['Bin']
        Completeness = float(row['Completeness'])
        Contamination = float(row['Contamination'])

        source = os.path.join(FastaDir,Bin)+'.fasta'
        if Completeness >= 50 and Contamination <= 10:
            target = os.path.join(Pass,Bin)+'.fasta'
        elif Completeness >= 50 and Contamination >10:
            target = os.path.join(Bandage, Bin) + '.fasta'
            if FastgDir is not None:
                shutil.copy(os.path.join(FastgDir,Bin)+'.fastg', os.path.join(Bandage,Bin)+'.fastg')
        else:
            target = os.path.join(Abandon, Bin) + '.fasta'

        shutil.copy(source, target)





'''
FastaQC2后，怎样输入fastg返回校正后的fasta文件需要用户自己手动操作
'''

