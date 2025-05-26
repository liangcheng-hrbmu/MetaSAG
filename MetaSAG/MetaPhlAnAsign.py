import math
import re
import sys
import subprocess
import os
import numpy as np
+2345
import pandas as pd
from collections import Counter
import warnings
import time
'''
脚本注释
1. 输入fastq文件，对fastq文件进行MetaPhlan4的注释 函数：MPAnno()
2. 对MetaPhlan4得到的bowtie2注释文件进行Cell_Align[KnownCell DoubleCell UnknownCell UnAsignCell]

'''

def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    return wrapper


######Step1 Metaphlan4注释fastq文件######
@timeit
def MPAnno(inputFastq, MPOut, SamName, env=None, DB=None):
    # 创建MPOut工作目录
    if not os.path.exists(MPOut):
        os.makedirs(MPOut, exist_ok=True)

    OutProfile = os.path.join(MPOut, SamName + '_profile.txt')
    OutBowtie = os.path.join(MPOut, SamName + '_bowtie2.bz2')
    OutVSC = os.path.join(MPOut, SamName + '_VSC.txt')

    if type(inputFastq) == str and inputFastq.endswith('.fastq'):
        ReadsEnd = 'Single'
        inputFastq = [inputFastq]
    elif type(inputFastq) == list:
        if len(inputFastq) == 2 and inputFastq[0].endswith('_R1.fastq') and inputFastq[1].endswith('_R2.fastq'):
            ReadsEnd = 'Pair'
            inputFastq1 = inputFastq[0]
            inputFastq2 = inputFastq[1]

        elif len(inputFastq) == 2 and inputFastq[0].endswith('_R2.fastq') and inputFastq[1].endswith('_R1.fastq'):
            ReadsEnd = 'Pair'
            inputFastq1 = inputFastq[1]
            inputFastq2 = inputFastq[0]

        elif len(inputFastq) == 1 and inputFastq[0].endswith('.fastq'):
            ReadsEnd = 'Single'
            inputFastq = inputFastq

        else:
            warnings.warn("传入的inputfastq不符合格式要求", UserWarning)

    else:
        warnings.warn("传入的inputfastq不符合格式要求", UserWarning)


    if ReadsEnd == 'Single':
        inputFastq = inputFastq[0]
    else:
        inputFastq = inputFastq1 + ',' + inputFastq2


    if DB is not None:
        command = f'metaphlan {inputFastq} --input_type fastq --profile_vsc --nproc 4 -o {OutProfile} --bowtie2out {OutBowtie} --vsc_out {OutVSC} --bowtie2db {DB}'
    else:
        command = f'metaphlan {inputFastq} --input_type fastq --profile_vsc --nproc 4 -o {OutProfile} --bowtie2out {OutBowtie} --vsc_out {OutVSC}'

    if env is not None:
        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)

    subprocess.call('bunzip2 ' + OutBowtie, shell=True)


######Step2 统计Bowtie2Anno文件中每个细胞 * 每个SGB的count情况######
# 得到四个文件Cell_SGB.txt Cell_Erro.txt Cell_SGB_Count.txt Cell_Erro_Count.txt

class MPBowtie():
    PYTHONDIR = os.path.dirname(os.path.abspath(__file__))

    def __init__(self,inputBowtie,outputDir):
        self.inputBowtie = inputBowtie
        self.outputDir = outputDir

    @timeit
    def CellSGBStatistic(self):

        InputBowtie = self.inputBowtie
        outputDir = self.outputDir

        # 创建outputDir
        if not os.path.exists(outputDir):
            os.makedirs(outputDir, exist_ok=True)

        Cell_SGB = os.path.join(outputDir, 'Cell_SGB.txt')
        Cell_Other = os.path.join(outputDir, 'Cell_Other.txt')
        Cell_SGB_Count = os.path.join(outputDir, 'Cell_SGB_Count.txt')
        Cell_Other_Count = os.path.join(outputDir, 'Cell_Other_Count.txt')
        AllCellSummary = os.path.join(outputDir, 'AllCellSummary.txt')

        # Step1 统计每个细胞 * 每个SGB的count
        with open(Cell_SGB, 'w') as output:
            output.write('Cell\tSGB\n')

        with open(Cell_Other, 'w') as output:
            output.write('Cell\tOther\n')

        SGB_freq_dict = {}
        Other_freq_dict = {}

        with open(InputBowtie) as input1:
            while True:
                line = input1.readline()
                if len(line) == 0:
                    break
                cell = line.split(':')[0]

                pattern = r"SGB\d{0,20}"
                SGB = re.search(pattern, line)
                if SGB is None:
                    with open(Cell_Other, 'a') as output:
                        anno = line.split('\t')[1]
                        anno = anno.strip()
                        output.write(cell + '\t' + anno + '\n')

                        if cell not in Other_freq_dict:
                            Other_freq_dict[cell] = Counter()

                        Other_freq_dict[cell][anno] += 1


                else:
                    with open(Cell_SGB, 'a') as output:
                        output.write(cell + '\t' + SGB.group() + '\n')
                        SGB = SGB.group()

                        if cell not in SGB_freq_dict:
                            SGB_freq_dict[cell] = Counter()

                        SGB_freq_dict[cell][SGB] += 1

        for cell, counter in SGB_freq_dict.items():
            for SGB, count in counter.items():
                with open(Cell_SGB_Count, 'a') as output:
                    output.write(str(cell) + '\t' + str(SGB) + '\t' + str(count) + '\n')

        for cell, counter in Other_freq_dict.items():
            for other, count in counter.items():
                with open(Cell_Other_Count, 'a') as output:
                    output.write(str(cell) + '\t' + str(other) + '\t' + str(count) + '\n')



        self.Cell_SGB_Count = pd.read_csv(Cell_SGB_Count, sep='\t', header=None)
        self.Cell_Other_Count = pd.read_csv(Cell_Other_Count, sep='\t', header=None)

    @timeit
    def CellAsign(self, BC_Count, Min_Reads=None ,total_annoreads = 200,SGB_rate = 0.8,Double_SGB_rate = 0.2, Double_All_rate = 0.8, Double_Annoreads = 200,Unknown_Allreads = 500, Unknown_Annoreads = 10, SGB_topCell=50):
        # 创建outputDir
        outputDir = self.outputDir
        if not os.path.exists(outputDir):
            os.makedirs(outputDir, exist_ok=True)

        AllCellSummary = os.path.join(outputDir, 'AllCellSummary.txt')
        KnownSGBFile = os.path.join(outputDir, 'KnownSGB.txt')
        DoubleCellFile = os.path.join(outputDir, 'DoubleCell.txt')
        UnknownSGBFile = os.path.join(outputDir, 'UnknownSGB.txt')
        UnAsignedSGBFile = os.path.join(outputDir, 'UnAsignedSGB.txt')
        KnownCellAssemFile = os.path.join(outputDir, 'KnownCellAssem_top50.txt')
        CellAnnoFile = os.path.join(outputDir, 'CellAnno.txt')

        # 1. 统计每个细胞——每个SGB——SGB.annoreads——total.annoreads——SGB.rate——allreads
        Cell_SGB_Count = self.Cell_SGB_Count
        #Cell_SGB_Count = pd.read_csv(Cell_SGB_Count, sep='\t', header=None)
        Cell_SGB_Count[3] = Cell_SGB_Count.groupby(Cell_SGB_Count.columns[0])[Cell_SGB_Count.columns[2]].transform('sum')
        Cell_SGB_Count[4] = Cell_SGB_Count[Cell_SGB_Count.columns[2]] / Cell_SGB_Count[Cell_SGB_Count.columns[3]]
        Cell_SGB_Count.columns = ['cell', 'SGB', 'SGB_annoreads', 'total_annoreads', 'SGB_rate']
        #bc_reads = pd.read_csv(BC_Count, sep='\t', header=0)
        bc_reads = BC_Count #BC_Count是BCFilterClass.py产生对象存储的值
        bc_reads.columns = ['cell', 'allreads']
        Cell_SGB_Count = pd.merge(Cell_SGB_Count, bc_reads[['cell', 'allreads']], on='cell', how='left')
        if Min_Reads is not None and Min_Reads > 0:
            Cell_SGB_Count = Cell_SGB_Count.loc[Cell_SGB_Count['allreads'] >= int(Min_Reads)]

        Cell_SGB_Count.to_csv(AllCellSummary, sep='\t', index=False)

        # 2. 分配已知物种的细胞 得到文件KnownSGB.file:cell   SGB total_annoreads  SGB_maxannorate cell_reads
        # annoreads >= 200 & SGB_maxannorate >= 0.8
        KnownCell = Cell_SGB_Count.loc[(Cell_SGB_Count['total_annoreads'] >= total_annoreads) & (Cell_SGB_Count['SGB_rate'] >= SGB_rate)]
        KnownCell.to_csv(KnownSGBFile, sep='\t', index=False)
        self.KnownCell = KnownCell

        # 3. 分配多物种液滴 得到文件DoubleCell.file:cell SGB total_annoreads SGB_annorate cell_reads
        # annoreads >=200 & Each_rate > 0.2 & All_rate >= 0.8
        DoubleCell = Cell_SGB_Count.loc[(Cell_SGB_Count['SGB_rate'] >= Double_SGB_rate) & (Cell_SGB_Count['total_annoreads'] >= Double_Annoreads)]
        DoubleCell['All_rate'] = DoubleCell.groupby('cell')['SGB_rate'].transform('sum')
        DoubleCell = DoubleCell.loc[DoubleCell['All_rate'] > Double_All_rate]
        count_temp = DoubleCell['cell'].value_counts()
        mask = DoubleCell['cell'].isin(count_temp[count_temp > 1].index)
        DoubleCell = DoubleCell.loc[mask]
        DoubleCell.to_csv(DoubleCellFile, sep='\t', index=False)
        self.DoubleCell = DoubleCell

        # 4. 分配未知物种液滴  得到文件UnknownCell.file: cell total_annoreads cell_reads
        # annoreads < 10 & allreads >1000 ()
        anno_zero = bc_reads.loc[~bc_reads['cell'].isin(Cell_SGB_Count['cell'])]
        UnknownCell = Cell_SGB_Count.loc[(Cell_SGB_Count['total_annoreads'] < Unknown_Annoreads) & (Cell_SGB_Count['allreads'] >= Unknown_Allreads)]

        if len(anno_zero) > 0:
            # print('type 1')
            anno_zero_bigcell = anno_zero['cell'].loc[anno_zero['allreads'] >= Unknown_Allreads]
            UnknownCell = UnknownCell['cell'].tolist()
            UnknownCell = list(set(UnknownCell))
            UnknownCell = UnknownCell + list(anno_zero_bigcell)
        else:
            # print('type 2')
            UnknownCell = UnknownCell['cell'].tolist()
            UnknownCell = list(set(UnknownCell))

        self.UnknownCell = UnknownCell

        with open(UnknownSGBFile, 'a') as output:
            for i in UnknownCell:
                output.write(str(i) + '\n')

        # 5. 分配除未知物种、已知物种、多胞液滴外的细胞及其最大SGB比例
        # (包括SGBannoreads==0的液滴、SGBannoreads!=0的液滴)

        AsignedCell = list(KnownCell['cell'].unique()) + list(DoubleCell['cell'].unique()) + list(set(UnknownCell))
        if Min_Reads is not None and Min_Reads > 0:
            minreads_bc = bc_reads.loc[bc_reads['allreads'] > int(Min_Reads)]
        else:
            minreads_bc = bc_reads

        CellMaxSGB = pd.DataFrame(columns=['cell', 'SGB', 'SGB_annoreads', 'total_annoreads', 'SGB_rate', 'allreads'])
        for Cell in list(Cell_SGB_Count['cell'].unique()):
            temp = Cell_SGB_Count.loc[Cell_SGB_Count['cell'] == Cell]
            max_SGB_row = temp.loc[temp['SGB_rate'].idxmax()]
            CellMaxSGB = CellMaxSGB._append(max_SGB_row, ignore_index=True)

        cell_temp2 = CellMaxSGB.loc[~CellMaxSGB['cell'].isin(AsignedCell)]
        anno_zero = bc_reads.loc[~bc_reads['cell'].isin(Cell_SGB_Count['cell'])]
        if len(anno_zero) > 0:
            cell_temp1 = anno_zero.loc[anno_zero['allreads'] < Unknown_Allreads]  # 没有任何SGB注释，且reads < 1000的细胞
            UnAsignedCell = pd.concat([cell_temp2, cell_temp1], ignore_index=True)
        else:
            UnAsignedCell = cell_temp2

        UnAsignedCell.to_csv(UnAsignedSGBFile, sep='\t', index=False)
        self.UnAsignedCell = UnAsignedCell
        # 6. 找出KnownSGB,找出每个SGB中前50个代表细胞
        KnownCellAssem = pd.DataFrame(
            columns=['cell', 'SGB', 'SGB_annoreads', 'total_annoreads', 'SGB_rate', 'allreads', 'total_annoreads_rank',
                     'allreads_rank', 'SGB_rate_rank', 'rank_sum', 'rank_rank'])

        for SGB in KnownCell['SGB'].unique():
            df_SGB = KnownCell[KnownCell['SGB'] == SGB]

            df_SGB['total_annoreads_rank'] = df_SGB['total_annoreads'].rank(method='dense', ascending=False)
            df_SGB['allreads_rank'] = df_SGB['allreads'].rank(method='dense', ascending=False)
            df_SGB['SGB_rate_rank'] = df_SGB['SGB_rate'].rank(method='dense', ascending=False)

            df_SGB['rank_sum'] = df_SGB['total_annoreads_rank'] + df_SGB['allreads_rank'] + df_SGB['SGB_rate_rank']
            df_SGB['rank_rank'] = df_SGB['rank_sum'].rank(method='dense', ascending=False)
            df_SGB_top50 = df_SGB.nlargest(SGB_topCell, 'rank_rank')

            KnownCellAssem = pd.concat([KnownCellAssem, df_SGB_top50], ignore_index=True)

        KnownCellAssem.to_csv(KnownCellAssemFile, sep='\t', index=False)
        self.KnownCellAssem = KnownCellAssem


        #将细胞分类整合为Cell_Anno
        cell_temp = self.KnownCell['cell'].tolist() + self.DoubleCell['cell'].unique().tolist() + self.UnknownCell + self.UnAsignedCell['cell'].tolist()
        type_temp = ['KnownCell'] * len(self.KnownCell['cell']) + ['DoubleCell'] * len(self.DoubleCell['cell'].unique()) + ['UnknownCell'] * len(self.UnknownCell) + ['UnAsignedCell'] * len(self.UnAsignedCell['cell'])
        sgb_temp = self.KnownCell['SGB'].tolist() + ['NoSGB'] * (len(cell_temp) - len(self.KnownCell['cell']))
        Cell_Anno = pd.DataFrame({'Cell':cell_temp,'Type':type_temp,'SGB':sgb_temp})
        Cell_Anno.to_csv(CellAnnoFile,index=False,sep='\t')
        self.Cell_Anno = Cell_Anno


    @timeit
    def HostPhage(self,Group='SGB',cumThresh=0.8 ,MinCellPhage=50):
        self.Cell_Other_Count.columns = ['Cell','Other','count']
        CellOther = self.Cell_Other_Count
        ResultDir = os.path.join(self.outputDir,'HostPhageResult')
        CellAnno = self.Cell_Anno


        # 创建目录
        if not os.path.exists(ResultDir):
            os.makedirs(ResultDir, exist_ok=True)

        frequency = CellOther['Other'].value_counts()
        frequency = pd.DataFrame(frequency)
        frequency_sorted = frequency.sort_values(by='count', ascending=False)
        frequency_sorted = pd.DataFrame(frequency_sorted)
        frequency_sorted = frequency_sorted.rename_axis('Other').reset_index()

        total = frequency['count'].sum()
        threshold = total * cumThresh

        # 找到达到80%阈值的行
        cumulative_sum = frequency_sorted['count'].cumsum()
        result_df = frequency_sorted[cumulative_sum <= threshold]
        MainPhage = result_df['Other'].tolist()

        # 结合CellErro,Cell_Anno构建cell phage group数据框
        CellPhageCluster = pd.merge(CellOther, CellAnno, on='Cell', how='left')
        PhageClusterCount = CellPhageCluster.groupby(['Other', Group]).size().reset_index(name='count')
        # PhageClusterCount.to_csv('./PhageClusterCount.txt',sep='\t',index=False)

        # table3 行Phage 列Cluster 元素Count
        unique_rows = PhageClusterCount['Other'].unique()
        unique_cols = PhageClusterCount[Group].unique()

        Table3 = np.zeros((len(unique_rows), len(unique_cols)))
        # 填充矩阵
        for index, row in PhageClusterCount.iterrows():
            row_idx = np.where(unique_rows == row['Other'])[0][0]
            col_idx = np.where(unique_cols == row[Group])[0][0]
            Table3[row_idx, col_idx] = row['count']

        Table3 = pd.DataFrame(Table3, index=unique_rows, columns=unique_cols)
        MainPhage_Table3 = Table3.loc[MainPhage, :]
        Table3.to_csv(os.path.join(ResultDir, 'Phage_Cluster_Count.txt'), sep='\t')

        self.Phage_Cluster_Count = Table3

        MainPhage_Table3.to_csv(os.path.join(ResultDir, 'MainPhage_Cluster_Count.txt'), sep='\t')

        self.MainPhage_Cluster_Count = MainPhage_Table3

        # table4 Phage Cluster Count_prop
        # 先统计每种Phage总的count有多少
        temp = PhageClusterCount.groupby('Other')['count'].sum().reset_index()
        Phage_count = dict(zip(temp['Other'], temp['count']))

        Table4 = pd.DataFrame(columns=Table3.columns)
        for index, row in Table3.iterrows():
            temp = row / Phage_count[index]
            Table4 = Table4._append(temp)
            # Table4.loc[len(Table4)] = temp

        MainPhage_Table4 = Table4.loc[MainPhage, :]

        Table4.to_csv(os.path.join(ResultDir, 'Phage_Cluster_CountProp.txt'), sep='\t')

        self.Phage_Cluster_CountProp = Table4

        MainPhage_Table4.to_csv(os.path.join(ResultDir, 'MainPhage_Cluster_CountProp.txt'), sep='\t')

        self.MainPhage_Cluster_CountProp = MainPhage_Table4

        # table1和table2 判断哪些细胞中含有噬菌体

        # 先统计cell phage count
        CellPhageCount = CellOther.groupby(['Other', 'Cell']).size().reset_index(name='count')
        CellPhageCount = pd.merge(CellPhageCount, CellAnno, on='Cell', how='left')
        # CellPhageCount.to_csv('./CellPhageCount.txt', sep='\t', index=False)

        # table1 行Phage 列Cluster 元素： 该Cluster下有多少细胞中含有Phage("含有"指该细胞中含有超过50条对应Phage的reads)
        temp = CellPhageCount[CellPhageCount['count'] >= MinCellPhage]

        unique_rows = temp['Other'].unique()
        unique_cols = temp[Group].unique()

        Table1 = np.zeros((len(unique_rows), len(unique_cols)))
        # 填充矩阵
        for index, row in temp.iterrows():
            row_idx = np.where(unique_rows == row['Other'])[0][0]
            col_idx = np.where(unique_cols == row[Group])[0][0]
            Table1[row_idx, col_idx] += 1

        Table1 = pd.DataFrame(Table1, index=unique_rows, columns=unique_cols)

        # table2 行Phage 列Cluster 元素：table1中细胞数，除以所在Cluster中细胞总数的比例
        Table2 = pd.DataFrame()

        temp = CellAnno.groupby(Group).size().reset_index()
        temp.columns = [Group, 'CellNum']
        Cluster_CellNum = dict(zip(temp[Group], temp['CellNum']))
        for colname in Table1.columns:
            Table2[colname] = Table1[colname] / Cluster_CellNum[colname]

        Table1.to_csv(os.path.join(ResultDir, 'Phage_Cluster_CellNum.txt'), sep='\t')

        self.Phage_Cluster_CellNum = Table1

        Table2.to_csv(os.path.join(ResultDir, 'Phage_Cluster_CellNumProp.txt'), sep='\t')

        self.Phage_Cluster_CellNumProp = Table2

    @timeit
    def DoubleCellKraken(self, CellBarn, env=None, KrakenDB=None, ReadsEnd='Single'):



        DoubleCellReport = os.path.join(self.outputDir, 'DoubleCellReport')
        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
        except ValueError as e:
            print(e)
            sys.exit()

        if not os.path.exists(DoubleCellReport):
            os.makedirs(DoubleCellReport, exist_ok=True)

        DoublePlotByPavian_R = os.path.join(MPBowtie.PYTHONDIR, 'DoublePlotByPavian.R')

        for dc in self.DoubleCell['cell'].unique():
            # kraken2注释每个细胞的fastq文件

            if ReadsEnd == 'Single':
                dcFastq = os.path.join(CellBarn, dc) + '.fastq'
                dcReport = os.path.join(DoubleCellReport, dc) + '_report'
                dcHtml = os.path.join(DoubleCellReport, dc) + '.html'
                if KrakenDB is not None:
                    command = 'kraken2 --db ' + KrakenDB + ' ' + dcFastq + ' --report ' + dcReport
                else:
                    command = 'kraken2 ' + dcFastq + ' --report ' + dcReport

                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

            elif ReadsEnd == 'Pair':
                dcFastq1 = os.path.join(CellBarn, dc) + '_R1.fastq'
                dcFastq2 = os.path.join(CellBarn, dc) + '_R2.fastq'
                dcReport = os.path.join(DoubleCellReport, dc) + '_report'
                dcHtml = os.path.join(DoubleCellReport, dc) + '.html'

                if KrakenDB is not None:
                    command = 'kraken2 --db ' + KrakenDB + ' ' + dcFastq1 + ' ' + dcFastq2 + ' --report ' + dcReport
                else:
                    command = 'kraken2 ' + dcFastq1 + ' ' + dcFastq2 + ' --report ' + dcReport

                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

            # pavian.R绘制每个细胞的kraken_report
            command = 'Rscript ' + DoublePlotByPavian_R + ' -i ' + dcReport + ' -o ' + dcHtml
            subprocess.call(command, shell=True)

    @timeit
    def CellAssem(self, CellBarn, env=None, ReadsEnd='Single'):
        # 创建目录
        BinFastq = os.path.join(self.outputDir, 'BinFastq')
        BinFasta = os.path.join(self.outputDir, 'BinFasta')
        BinFastg = os.path.join(self.outputDir, 'BinFastg')
        DirCreat = [BinFastq, BinFasta, BinFastg]
        for dir in DirCreat:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        #KnownCellAssem = pd.read_csv(KnownCellAssem, sep='\t', header=0)
        KnownCellAssem = self.KnownCellAssem
        # 将细胞复制到相应SGB文件夹中
        for index, row in KnownCellAssem.iterrows():
            cell = row['cell']

            try:
                if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                    raise ValueError('ReadsEnd must be Single or Pair!')
            except ValueError as e:
                print(e)
                sys.exit()

            if ReadsEnd == 'Single':
                cellFile = os.path.join(CellBarn, row['cell']) + '.fastq'
                SGB = row['SGB']
                SGBFastqFile = os.path.join(BinFastq, SGB) + '.fastq'
                subprocess.call('cat ' + cellFile + ' >> ' + SGBFastqFile, shell=True)

            elif ReadsEnd == 'Pair':
                cellFile1 = os.path.join(CellBarn, row['cell']) + '_R1.fastq'
                cellFile2 = os.path.join(CellBarn, row['cell']) + '_R2.fastq'
                SGB = row['SGB']
                SGBFastqFile1 = os.path.join(BinFastq, SGB) + '_R1.fastq'
                SGBFastqFile2 = os.path.join(BinFastq, SGB) + '_R2.fastq'
                subprocess.call('cat ' + cellFile1 + ' >> ' + SGBFastqFile1, shell=True)
                subprocess.call('cat ' + cellFile2 + ' >> ' + SGBFastqFile2, shell=True)

        if ReadsEnd == 'Single':
            for SGB in KnownCellAssem['SGB'].unique():
                SGBFastqFile = os.path.join(BinFastq, SGB) + '.fastq'
                SGBFastaFile = os.path.join(BinFasta, SGB) + '.fasta'
                SGBFastgFile = os.path.join(BinFastg, SGB) + '.fastg'

                command = 'spades.py --sc --pe1-s ' + SGBFastqFile + ' -o  ' + BinFasta + '/' + SGB + '_sc'

                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

                subprocess.call('cp ' + BinFasta + '/' + SGB + '_sc/contigs.fasta ' + SGBFastaFile, shell=True)
                subprocess.call('cp ' + BinFasta + '/' + SGB + '_sc/assembly_graph.fastg ' + SGBFastgFile, shell=True)

        elif ReadsEnd == 'Pair':
            for SGB in KnownCellAssem['SGB'].unique():
                SGBFastqFile1 = os.path.join(BinFastq, SGB) + '_R1.fastq'
                SGBFastqFile2 = os.path.join(BinFastq, SGB) + '_R2.fastq'
                SGBFastaFile = os.path.join(BinFasta, SGB) + '.fasta'
                SGBFastgFile = os.path.join(BinFastg, SGB) + '.fastg'

                command = 'spades.py --sc --pe1-1 ' + SGBFastqFile1 + ' --pe1-2 ' + SGBFastqFile2 + ' -o ' + BinFasta + '/' + SGB + '_sc'

                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

                subprocess.call('cp ' + BinFasta + '/' + SGB + '_sc/contigs.fasta ' + SGBFastaFile, shell=True)
                subprocess.call('cp ' + BinFasta + '/' + SGB + '_sc/assembly_graph.fastg ' + SGBFastgFile, shell=True)


'''
from zszshhh import MetaPhlanAsignClass
import pandas as pd
x=MetaPhlanAsignClass.MPBowtie('D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\MetaPhlanAsign\\input\\toy.bowtie2','D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\MetaPhlanAsign\\result')
x.inputBowtie
x.outputDir
x.CellSGBStatistic()
x.Cell_SGB_Count
x.Cell_Erro_Count
BC_Count=pd.read_csv('D:\\Tools\\PyCharm 2023.3.4\\MicrobePlusTest\\TestPyPI\\zszshhh\\testData\\MetaPhlanAsign\\input\\bc_read',sep='\t',header=None)
x.CellAsign(BC_Count)
x.KnownCell
x.DoubleCell
x.UnknownCell
x.UnAsignedCell
x.KnownCellAssem
x.HostPhage()
'''


