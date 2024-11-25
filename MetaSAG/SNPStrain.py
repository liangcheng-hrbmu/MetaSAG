import shutil

import pandas as pd
import os
import subprocess
import json
import numpy as np
import umap
from sklearn import preprocessing
import time
import pandas as pd
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree, fcluster
import matplotlib
import matplotlib.pyplot as plt
import os
import random
import matplotlib.lines as mlines
import matplotlib as mpl
from pylab import *
import json
from ast import literal_eval
from matplotlib.cm import ScalarMappable
import scipy
import scipy.stats
import scipy.signal
from matplotlib.ticker import MaxNLocator
import subprocess
import sys
import pickle

'''
脚本注释：
分菌株SNP脚本提供两个策略：
1. 每次只对一个物种的细胞箱进行分菌株。该细胞箱中应该有每个细胞的fastq文件，以及该箱的组装基因组文件
2. 一次对所有箱的细胞进行分菌株。要求统一提供组装基因组文件以及存放每个细胞fastq文件的CellBarn文件夹路径。[提供可选参数，表示对所有细胞中哪些箱进行分菌株的操作]
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












class SingleBin():

    def __init__(self,BinDir, ResultDir, SpeciesName):
        self.BinDir = BinDir
        self.ResultDir = ResultDir
        self.SpeciesName = SpeciesName


    strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                         "#8c564b",
                         "#e377c2"]
    umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                          u'#7f7f7f',
                          u'#bcbd22', u'#17becf']

    def strain_split1(self,work_dirt_temp='', file_name='', cutoff_cell_ratio=0.05, cutoff_type_count=2,
                      cutoff_type_ratio=0.01,
                      cutoff_SNP_ratio=0.01):
        barcode_list = []
        contig_list = []
        location_on_contig_list = []
        reference_list = []
        alter_list = []
        sequence_list = []

        total_line = 0
        with open(os.path.join(work_dirt_temp, file_name)) as finput:
            while True:
                line_read = finput.readline()
                if len(line_read) == 0:
                    break
                total_line += 1
                if total_line == 4:
                    parsed_bam = line_read.split()
                    for i in range(len(parsed_bam)):
                        barcode_with_bam = parsed_bam[i]
                        if '.bam' == barcode_with_bam[-4:]:
                            barcode_list += [barcode_with_bam[:-4]]
                        cutoff_cell_count = int(cutoff_cell_ratio * len(barcode_list))
                        cutoff_type_count = max(cutoff_type_count, int(cutoff_type_ratio * len(barcode_list)))

                if line_read[0] != '#':
                    line_split_temp = line_read.split()
                    contig_list += [int(line_split_temp[0].split('_')[1])]
                    location_on_contig_list += [int(line_split_temp[1])]
                    reference_list += [line_split_temp[3]]
                    alter_list += [line_split_temp[4]]

                    mapping_count_list_temp = []
                    for i in range(9, len(line_split_temp)):
                        mapping_result_temp = line_split_temp[i].split(':')[-1]
                        reference_count_temp = int(mapping_result_temp.split(',')[0])
                        alter_count_temp = int(mapping_result_temp.split(',')[1])
                        total_count_temp = reference_count_temp + alter_count_temp
                        if total_count_temp < 2:
                            mapping_count_list_temp += [0]
                        elif max(reference_count_temp, alter_count_temp) / float(total_count_temp) < 0.99:
                            mapping_count_list_temp += [0]
                        else:
                            if max(reference_count_temp, alter_count_temp) == reference_count_temp:
                                mapping_count_list_temp += [1]
                            elif max(reference_count_temp, alter_count_temp) == alter_count_temp:
                                mapping_count_list_temp += [-1]
                            else:
                                print('someting wrong')
                    sequence_list += [mapping_count_list_temp]

        #print(len(sequence_list))
        #print(sequence_list[0])

        snp_pd = pd.DataFrame({'1': sequence_list[0]})
        for i in range(len(sequence_list)):
            snp_pd[str(i + 1)] = sequence_list[i]
            # snp_pd=snp_pd.copy()
        snp_pd = snp_pd.T
        #print(snp_pd)
        #print(barcode_list)
        snp_pd.columns = barcode_list

        #print(snp_pd.shape)
        #print(len(barcode_list))
        #print(barcode_list)

        snp_pd['contig_list'] = contig_list
        snp_pd['location_on_contig_list'] = location_on_contig_list
        snp_pd['reference_list'] = reference_list
        snp_pd['alter_list'] = alter_list
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp.csv')

        to_drop_row_list = []
        for i in range(len(snp_pd)):
            if sum(snp_pd.iloc[i][:-4] != 0) < cutoff_cell_count:
                to_drop_row_list.append(i)
                continue
            if sum(snp_pd.iloc[i][:-4] == -1) < cutoff_type_count or sum(snp_pd.iloc[i][:-4] == 1) < cutoff_type_count:
                to_drop_row_list.append(i)
        print('File name is {},total SNPs before drop is {},SNPs to drop is {}.'.format(file_name, len(snp_pd),
                                                                                        len(to_drop_row_list)))

        self.DropSNPNum = len(to_drop_row_list)
        self.AllSNPNum = len(snp_pd)

        snp_pd = snp_pd.drop(to_drop_row_list, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered.csv')
        cutoff_SNP_count = max(cutoff_SNP_ratio * len(snp_pd), 10)
        to_drop_column_list = []
        barcode_list = list(snp_pd.columns.values[:-4])
        for i in range(snp_pd.shape[1] - 4):
            if sum(snp_pd[barcode_list[i]] != 0) < cutoff_SNP_count:
                to_drop_column_list += [barcode_list[i]]

        print("Total cells before drop is {},cell to drop is {}.".format(snp_pd.shape[1] - 4, len(to_drop_column_list)))
        self.AllCellNum = snp_pd.shape[1] - 4
        self.DropCellNum = len(to_drop_column_list)

        snp_pd = snp_pd.drop(to_drop_column_list, axis=1)
        snp_pd.to_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered_bc_dropped.csv', index=False)

        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name) + '_snp_filtered_bc_dropped.csv')
        total_ratio_of_snps_list = []
        for i in range(len(snp_pd)):
            total_coverage_temp = sum(list(snp_pd.iloc[i][:-4] != 0))
            if total_coverage_temp == 0:
                print(i)
                total_ratio_of_snps_list += [1]
                continue
            total_ratio_of_snps_list += [sum(list(snp_pd.iloc[i][:-4] == 1)) / float(total_coverage_temp)]
            plt.figure()
            count_temp = plt.hist(np.array(total_ratio_of_snps_list), bins=np.arange(0, 1.01, 0.02))
            plt.tick_params(axis='both', direction='in', which='both', top=True, right=True)
            plt.xlabel('Ratio of cells matching strain 1')
            plt.ylabel('Number of alleles')
            plt.savefig(os.path.join(work_dirt_temp, file_name) + '_SNP.svg', format='svg', bbox_inches='tight')
            plt.clf()
            return (total_ratio_of_snps_list)

    def check_strain1(self,file_name, dirt_temp, output_dirt, tree_cut_threshold=2.7):
        strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                             "#8c564b", "#e377c2"]
        umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
        default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                              u'#7f7f7f', u'#bcbd22', u'#17becf']

        # work_dirt_temp=dirt_base+'fastq/'+bin+'/bamResult/variants/'
        work_dirt_temp = dirt_temp
        snp_pd = pd.read_csv(os.path.join(work_dirt_temp, file_name))
        #print(snp_pd.shape)  # 行是SNP，列是CELL
        snp_pd_temp = snp_pd
        snp_pd_temp = snp_pd_temp.drop(['contig_list', 'location_on_contig_list', 'reference_list', 'alter_list'],
                                       axis=1)
        snp_pd_temp = snp_pd_temp.T
        # print(snp_pd_temp)
        X = preprocessing.normalize(snp_pd_temp)
        linked = linkage(X, method='ward')
        # print(linked)
        matplotlib.rcParams.update({'font.size': 50})
        fig, ax = plt.subplots(1, 1, figsize=(3, 9))

        plt.yticks([0, 50, 100, 150, 200])
        plt.tick_params(axis='both', direction="in", which='both', right=True, length=15, width=3)
        matplotlib.rcParams['lines.linewidth'] = 0.8

        # strain_color_list_BV = [strain_color_list[3], strain_color_list[2], strain_color_list[1],strain_color_list[0]]###############################################################
        strain_color_list_BV = [strain_color_list[2], strain_color_list[0]]
        hierarchy.set_link_color_palette(strain_color_list_BV)

        # tree_cut_threshold=1.8
        dendrogram_plot = dendrogram(linked,
                                     orientation='left',
                                     distance_sort='ascending',
                                     leaf_font_size=8,
                                     no_labels=True,
                                     color_threshold=tree_cut_threshold,
                                     above_threshold_color='black',
                                     show_leaf_counts=False)

        plt.xticks([])
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_tree.svg',
                    format='svg',
                    bbox_inches='tight')
        matplotlib.rcParams.update({'font.size': 12})

        plt.figure()
        X = snp_pd_temp
        X = preprocessing.normalize(snp_pd_temp)
        reducer = umap.UMAP(random_state=42)
        embedding = reducer.fit_transform(X)
        #print(embedding.shape)
        # print(embedding)
        plt.scatter(embedding[:, 0], embedding[:, 1], s=3)
        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_umap.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()

        index_int_list = []

        for i in dendrogram_plot['ivl']:
            index_int_list.append(int(i))
        SNP_similarity_df = pd.DataFrame(-1.0, index=snp_pd.columns[index_int_list],
                                         columns=snp_pd.columns[index_int_list])
        for i in range(len(SNP_similarity_df)):
            column_name_a = snp_pd.columns[index_int_list[i]]
            for j in range(len(SNP_similarity_df)):
                if i == j:
                    SNP_similarity_df.iloc[i][j] = 1.0
                    continue
                column_name_b = snp_pd.columns[index_int_list[j]]
                check_snp_temp = snp_pd[[column_name_a, column_name_b]]
                check_snp_temp = check_snp_temp.loc[check_snp_temp[column_name_b] != 0]
                check_snp_temp = check_snp_temp.loc[check_snp_temp[column_name_a] != 0]

                if len(check_snp_temp) == 0:
                    continue
                SNP_similarity_df.iloc[i][j] = sum(
                    check_snp_temp[column_name_b] == check_snp_temp[column_name_a]) / float(
                    len(check_snp_temp))

        SNP_similarity_df.to_csv(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_df.csv',
                                 index=False)
        fig, ax = plt.subplots()
        data_color = []
        fig.set_size_inches(6, 6, forward=True)
        # my_cmap=plt.cm.get_cmap(matplotlib.cm.cmap_d.keys()[58])
        # my_cmap = mpl.colors.ListedColormap(['w', default_color_list[4]])
        my_cmap = 'RdBu'
        maximum_similarity = 1

        # print(matplotlib.cm.cmap_d.keys()[58])
        SNP_similarity_df = pd.read_csv(
            os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity_df.csv')
        ax.imshow(SNP_similarity_df, cmap=my_cmap, vmax=maximum_similarity, vmin=0)
        ax.set_xlabel("Single cells", fontsize=15)
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0, maximum_similarity))
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='vertical', fraction=0.05, pad=0.05, aspect=15,
                            ticks=[0, maximum_similarity])
        cbar.set_label('Similarity of SNPs between cells', rotation=90, labelpad=-1, fontsize=15)
        plt.savefig(os.path.join(output_dirt, file_name.split('_')[0]) + 'sp_call_SNP_similarity.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()
        return (snp_pd, linked, embedding)

    # 要求SpeciesName中不能有'_'
    @timeit
    def SingleBinPrepare(self, bcftools=None, snap_aligner=None, env=None, ReadsEnd='Single'):

        BinDir = self.BinDir
        ResultDir = self.ResultDir
        SpeciesName = self.SpeciesName

        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
        except ValueError as e:
            print(e)
            sys.exit()

        # 创建目录
        ResultBam = os.path.join(ResultDir, 'ResultBam')
        ResultVCF = os.path.join(ResultDir, 'ResultVCF')
        ResultPic = os.path.join(ResultDir, 'ResultPic')

        self.PrepareBam = ResultBam
        self.PrepareVCF = ResultVCF
        self.PreparePic = ResultPic

        DirCreate = [ResultDir, ResultBam, ResultVCF, ResultPic]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        if bcftools is None:
            bcftools = 'bcftools'

        if snap_aligner is None:
            snap_aligner = 'snap-aligner'

        FastaFile = os.path.join(BinDir, SpeciesName) + '.fasta'
        FastaIndex = os.path.join(ResultBam, SpeciesName) + '_index'
        ###### 对Bin的fasta构建索引[SNAP] ######
        command = snap_aligner + ' index ' + FastaFile + ' ' + FastaIndex
        subprocess.run(['conda', 'run', '-n', 'base', 'bash', '-c', command])

        ###### 比对得到bam vcf文件[SNAP] ######
        FastqList = []
        if ReadsEnd == 'Single':
            for filename in os.listdir(BinDir):
                if filename.endswith('.fastq'):
                    FastqList.append(filename[:-6])

            for cell in FastqList:
                CellFastq = os.path.join(BinDir, cell) + '.fastq'
                CellBam = os.path.join(ResultBam, cell) + '.bam'
                command = snap_aligner + ' single ' + FastaIndex + ' ' + CellFastq + ' -I -R ID:FASTQ\\tPL:Illumina\\t:PU:pu\\tLB:lb\\tSM:' + cell + ' -t 8 -so -o ' + CellBam
                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

        elif ReadsEnd == 'Pair':
            for filename in os.listdir(BinDir):
                if filename.endswith('_R1.fastq'):
                    FastqList.append(filename[:-9])

            for cell in FastqList:
                CellFastq_R1 = os.path.join(BinDir, cell) + '_R1.fastq'
                CellFastq_R2 = os.path.join(BinDir, cell) + '_R2.fastq'
                CellBam = os.path.join(ResultBam, cell) + '.bam'
                command = snap_aligner + ' paired ' + FastaIndex + ' ' + CellFastq_R1 + ' ' + CellFastq_R2 + ' -I -R ID:FASTQ\\tPL:Illumina\\t:PU:pu\\tLB:lb\\tSM:' + cell + ' -t 8 -so -o ' + CellBam
                if env is not None:
                    subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

        ###### 使用bcftools获取vcf文件[bcftools] ######
        command = bcftools + ' mpileup -f ' + FastaFile + ' -a FORMAT/AD ' + ResultBam + '/*.bam > ' + ResultVCF + '/raw_calls.bcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + ' call --ploidy 1 -v -m ' + ResultVCF + '/raw_calls.bcf > ' + ResultVCF + '/calls.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + ' view -v snps ' + ResultVCF + '/calls.vcf > ' + ResultVCF + '/calls_snps.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        command = bcftools + " filter -o " + ResultVCF + '/' + SpeciesName + "_calls_snps_qual30.vcf -i 'QUAL>30' " + ResultVCF + '/calls_snps.vcf'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

        ######
        self.strain_split1(ResultVCF, SpeciesName + '_calls_snps_qual30.vcf')
        snp_pd, linked, embedding = self.check_strain1(SpeciesName + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv',
                                                  ResultVCF, ResultPic)
        plt.clf()

        with open(ResultVCF + '/' + 'snp_pd.pkl', 'wb') as f:
            pickle.dump(snp_pd, f)

        with open(ResultVCF + '/' + 'linked.pkl', 'wb') as f:
            pickle.dump(linked, f)

        with open(ResultVCF + '/' + 'embedding.pkl', 'wb') as f:
            pickle.dump(embedding, f)

        self.snp_pd = snp_pd
        self.linked = linked
        self.embedding = embedding

    @timeit
    def SingleBinSplit(self,ClusterNum):
        ClusterNum = int(ClusterNum)
        Species = self.SpeciesName
        PrepareDir = self.ResultDir

        PrepareBam = self.PrepareBam
        PrepareVCF = self.PrepareVCF
        PreparePic = self.PreparePic


        #with open(PrepareVCF + '/linked.pkl', 'rb') as f:
            #linked = pickle.load(f)

        #with open(PrepareVCF + '/snp_pd.pkl', 'rb') as f:
            #snp_pd = pickle.load(f)

        #with open(PrepareVCF + '/embedding.pkl', 'rb') as f:
            #embedding = pickle.load(f)



        snp_pd = self.snp_pd
        linked = self.linked
        embedding = self.embedding

        cutree_contigs = cut_tree(linked, n_clusters=ClusterNum)
        umap_compare = pd.DataFrame({'Barcode': snp_pd.columns[:-4]})
        umap_compare['X'] = embedding[:, 0]
        umap_compare['Y'] = embedding[:, 1]

        cluster_hier = []
        for i in range(len(cutree_contigs)):
            cluster_hier.append(cutree_contigs[i][0])

        umap_compare['cluster_hier'] = cluster_hier

        for j in range(ClusterNum):
            umap_compare_bin = umap_compare.loc[umap_compare['cluster_hier'] == j]
            plt.scatter(umap_compare_bin['X'], umap_compare_bin['Y'], s=10, c=SingleBin.strain_color_list[j],
                        marker=SingleBin.umap_plot_shape[j])

        plt.xticks([])
        plt.yticks([])
        plt.ylabel('UMAP dimention 2')
        plt.xlabel('UMAP dimention 1')
        plt.savefig(os.path.join(PreparePic, Species) + '_sp_call_SNP_umap_final.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()

        ###################################
        ###get strain share SNP genotype###
        ###################################

        # file = os.path.join(PrepareVCF,Species) + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv'

        for i in range(ClusterNum):
            index_of_this_cluster = []
            for j in range(len(cutree_contigs)):
                if cutree_contigs[j] == i:
                    index_of_this_cluster += [j]
                snp_pd['hierarchy_cluster_' + str(i) + 'pos'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == 1).sum(axis=1)
                snp_pd['hierarchy_cluster_' + str(i) + 'neg'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == -1).sum(axis=1)

        strain_X_genotype = [[] for i in range(ClusterNum)]  # 创建长度为ClusterNum的空二维数组
        for i in range(len(snp_pd)):
            for j in range(ClusterNum):
                # generate genotype of each strain
                # at each SNP, if both have less than 2 cells, assign it 0, means not sure
                # if the max genotype is less than 90% total, assign 0
                # otherwise assign the major genotype
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] <= 1 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][
                            i] <= 1:
                    strain_X_genotype[j].append(0)
                    continue
                elif snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] != 0 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i] != 0:
                    if max(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                           snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]) / float(
                        min(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                            snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i])) < 9:
                        strain_X_genotype[j].append(0)
                        continue
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] > snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]:
                    strain_X_genotype[j].append(1)
                else:
                    strain_X_genotype[j].append(-1)

        for j in range(ClusterNum):
            snp_pd['strain_' + str(j) + '_genotype'] = strain_X_genotype[j]

        rows_to_drop = []
        for i in range(len(snp_pd)):
            genotype_SNP_temp = []
            for j in range(ClusterNum):
                genotype_SNP_temp.append(snp_pd['strain_' + str(j) + '_genotype'][i])
            if not (1 in genotype_SNP_temp and -1 in genotype_SNP_temp):
                rows_to_drop.append(i)

        print("SNPs before drop is {}, SNPs to drop is {}.".format(len(snp_pd), len(rows_to_drop)))

        snp_pd = snp_pd.drop(rows_to_drop, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)

        # this is to generate strain list
        list_of_ratio = [[] for i in range(ClusterNum)]
        barcode_list_strain = [[] for i in range(ClusterNum + 1)]

        for i in range(snp_pd.shape[1]):
            column_name_temp = snp_pd.columns[i]
            if column_name_temp.startswith('strain') or column_name_temp.startswith('hierarchy_cluster'):
                continue

            if 'list' not in column_name_temp:
                y = column_name_temp.split('/')[-1]
            else:
                continue

            # if not column_name_temp in snp_pd.columns:
            # continue
            temp_df_check = snp_pd.loc[snp_pd[column_name_temp] != 0]
            for j in range(ClusterNum):
                temp_df_check_strain = temp_df_check.loc[temp_df_check['strain_' + str(j) + '_genotype'] != 0]
                total_count_same = sum(temp_df_check_strain[column_name_temp] == \
                                       temp_df_check_strain['strain_' + str(j) + '_genotype'])
                list_of_ratio[j].append(total_count_same / float(len(temp_df_check_strain[column_name_temp])))

            list_of_ratio_check = [list_of_ratio[j][-1] for j in range(ClusterNum)]
            print(list_of_ratio_check)

            if max(list_of_ratio_check) < 0.95:
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            elif sum(np.array(list_of_ratio_check) > 0.95) > 1:
                print('Not sure,more then two strains match', column_name_temp)
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            else:
                # barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(column_name_temp)
                barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(y)

        with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'w') as output:
            output.write('Cluster\tCell\n')

        for cluster in range(len(barcode_list_strain)):
            # print('look here is barcode_list_strain:', i)
            # print(barcode_list_strain[i])
            for cell in barcode_list_strain[cluster]:
                with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                    output.write(str(cluster) + '\t' + str(cell) + '\n')

'''
CellAnno文件
 
Cluster Cell
SGB10068        Sam1025_10158
SGB10068        Sam1025_10251
SGB10068        Sam1025_10333
'''






class AllBin():
    def __init__(self,FastaDir, FastqDir, CellAnno, ResultDir):
        self.FastaDir = FastaDir
        self.FastqDir = FastqDir
        self.CellAnno = CellAnno
        self.ResultDir = ResultDir

    strain_color_list = ["#33bbff", "#0050e6", "#009999", "#777777", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                         "#8c564b",
                         "#e377c2"]
    umap_plot_shape = ['D', 'o', 's', 'v', 'd', 'p', '*', 'H', '+', 'x']
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2',
                          u'#7f7f7f',
                          u'#bcbd22', u'#17becf']
    
    @timeit
    def AllBinPrepare(self, bcftools=None, snap_aligner=None, env=None, ReadsEnd='Single'):

        FastaDir = self.FastaDir
        FastqDir = self.FastqDir
        CellAnno = self.CellAnno
        ResultDir = self.ResultDir

        ######  将所有物种细胞整理成SingleBinPrepare所需文件的目录格式  ######
        # ResultDir：ClusterxxxInput ClusterxxxResult（要求CellAnno文件中Cluster的名字与fasta文件名一致）
        try:
            if not (ReadsEnd == 'Single' or ReadsEnd == 'Pair'):
                raise ValueError('ReadsEnd must be Single or Pair!')
        except ValueError as e:
            print(e)
            sys.exit()

        CellAnno = pd.read_csv(CellAnno, sep='\t', header=0)
        if bcftools is None:
            bcftools = 'bcftools'

        if snap_aligner is None:
            snap_aligner = 'snap-aligner'

        # 创建所有ClusterInput文件夹
        DirCreate = []
        for cluster in CellAnno['Cluster'].unique():
            dir1 = os.path.join(ResultDir, str(cluster)) + 'Input'
            dir2 = os.path.join(ResultDir, str(cluster)) + 'Result'
            DirCreate.append(dir1)
            DirCreate.append(dir2)

        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        # 移动细胞的fastq以及箱的fasta文件到相应的ClusterxxxInput文件夹中
        if ReadsEnd == 'Pair':
            for index, row in CellAnno.iterrows():
                Cluster = str(row['Cluster'])
                Cell = str(row['Cell'])

                sourceFastq1 = os.path.join(FastqDir, Cell) + '_R1.fastq'
                sourceFastq2 = os.path.join(FastqDir, Cell) + '_R2.fastq'
                targetFastq1 = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '_R1.fastq'
                targetFastq2 = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '_R2.fastq'

                shutil.copy(sourceFastq1, targetFastq1)
                shutil.copy(sourceFastq2, targetFastq2)

        elif ReadsEnd == 'Single':
            for index, row in CellAnno.iterrows():
                Cluster = str(row['Cluster'])
                Cell = str(row['Cell'])

                sourceFastq = os.path.join(FastqDir, Cell) + '.fastq'
                targetFastq = os.path.join(ResultDir, str(Cluster) + 'Input', Cell) + '.fastq'

                shutil.copy(sourceFastq, targetFastq)

        # 同时为下面的AllBinSplit()函数准备BinClusterAnno文件
        # 需要输入想分菌株的ClusterxxxResult路径以及该箱判断分为几个簇。
        # BinClusterAnno.txt (SpeciesName必须与箱下的fasta文件名一致)
        # BinPrepareDir    ClusterNum    SpeciesName
        # /data_alluser/../SGB4933Result 4(分簇个数不能超过10)   SGB4933

        self.BinClusterAnno = os.path.join(ResultDir,'BinClusterAnno.txt')
        with open(os.path.join(ResultDir,'BinClusterAnno.txt'),'w') as f:
            f.write('BinPrepareDir\tSpeciesName\tClusterNum\n')


        for Cluster in CellAnno['Cluster'].unique():
            sourceFasta = os.path.join(FastaDir, str(Cluster)) + '.fasta'
            targetFasta = os.path.join(ResultDir, str(Cluster) + 'Input', str(Cluster)) + '.fasta'
            shutil.copy(sourceFasta, targetFasta)
            BinDir = os.path.join(ResultDir, str(Cluster)) + 'Input'
            ResultDirSingle = os.path.join(ResultDir, str(Cluster)) + 'Result'
            with open(os.path.join(ResultDir,'BinClusterAnno.txt'),'a') as f:
                f.write(ResultDirSingle + '\t' + str(Cluster) +'\n')

            sbp = SingleBin(BinDir,ResultDirSingle,str(Cluster))
            #SingleBinPrepare(BinDir, ResultDirSingle, str(Cluster), bcftools=bcftools, snap_aligner=snap_aligner, env=env, ReadsEnd=ReadsEnd)
            sbp.SingleBinPrepare(bcftools=bcftools, snap_aligner=snap_aligner, env=env, ReadsEnd=ReadsEnd)


    def SingleBinSplit(self,PrepareDir, ClusterNum, Species):
        ClusterNum = int(ClusterNum)

        # PrepareBam = os.path.join(PrepareDir, 'ResultBam')
        PrepareVCF = os.path.join(PrepareDir, 'ResultVCF')
        PreparePic = os.path.join(PrepareDir, 'ResultPic')

        with open(PrepareVCF + '/linked.pkl', 'rb') as f:
            linked = pickle.load(f)

        with open(PrepareVCF + '/snp_pd.pkl', 'rb') as f:
            snp_pd = pickle.load(f)

        with open(PrepareVCF + '/embedding.pkl', 'rb') as f:
            embedding = pickle.load(f)

        cutree_contigs = cut_tree(linked, n_clusters=ClusterNum)
        umap_compare = pd.DataFrame({'Barcode': snp_pd.columns[:-4]})
        umap_compare['X'] = embedding[:, 0]
        umap_compare['Y'] = embedding[:, 1]

        cluster_hier = []
        for i in range(len(cutree_contigs)):
            cluster_hier.append(cutree_contigs[i][0])

        umap_compare['cluster_hier'] = cluster_hier

        for j in range(ClusterNum):
            umap_compare_bin = umap_compare.loc[umap_compare['cluster_hier'] == j]
            plt.scatter(umap_compare_bin['X'], umap_compare_bin['Y'], s=10, c=AllBin.strain_color_list[j],
                        marker=AllBin.umap_plot_shape[j])

        plt.xticks([])
        plt.yticks([])
        plt.ylabel('UMAP dimention 2')
        plt.xlabel('UMAP dimention 1')
        plt.savefig(os.path.join(PreparePic, Species) + '_sp_call_SNP_umap_final.svg', format='svg',
                    bbox_inches='tight')
        plt.clf()

        ###################################
        ###get strain share SNP genotype###
        ###################################

        # file = os.path.join(PrepareVCF,Species) + '_calls_snps_qual30.vcf_snp_filtered_bc_dropped.csv'

        for i in range(ClusterNum):
            index_of_this_cluster = []
            for j in range(len(cutree_contigs)):
                if cutree_contigs[j] == i:
                    index_of_this_cluster += [j]
                snp_pd['hierarchy_cluster_' + str(i) + 'pos'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == 1).sum(axis=1)
                snp_pd['hierarchy_cluster_' + str(i) + 'neg'] = (
                        snp_pd[snp_pd.columns.values[index_of_this_cluster]] == -1).sum(axis=1)

        strain_X_genotype = [[] for i in range(ClusterNum)]  # 创建长度为ClusterNum的空二维数组
        for i in range(len(snp_pd)):
            for j in range(ClusterNum):
                # generate genotype of each strain
                # at each SNP, if both have less than 2 cells, assign it 0, means not sure
                # if the max genotype is less than 90% total, assign 0
                # otherwise assign the major genotype
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] <= 1 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][
                            i] <= 1:
                    strain_X_genotype[j].append(0)
                    continue
                elif snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] != 0 and \
                        snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i] != 0:
                    if max(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                           snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]) / float(
                        min(snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i],
                            snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i])) < 9:
                        strain_X_genotype[j].append(0)
                        continue
                if snp_pd['hierarchy_cluster_' + str(j) + 'pos'][i] > snp_pd['hierarchy_cluster_' + str(j) + 'neg'][i]:
                    strain_X_genotype[j].append(1)
                else:
                    strain_X_genotype[j].append(-1)

        for j in range(ClusterNum):
            snp_pd['strain_' + str(j) + '_genotype'] = strain_X_genotype[j]

        rows_to_drop = []
        for i in range(len(snp_pd)):
            genotype_SNP_temp = []
            for j in range(ClusterNum):
                genotype_SNP_temp.append(snp_pd['strain_' + str(j) + '_genotype'][i])
            if not (1 in genotype_SNP_temp and -1 in genotype_SNP_temp):
                rows_to_drop.append(i)

        print("SNPs before drop is {}, SNPs to drop is {}.".format(len(snp_pd), len(rows_to_drop)))
        snp_pd = snp_pd.drop(rows_to_drop, axis=0)
        snp_pd = snp_pd.reset_index(drop=True)

        # this is to generate strain list
        list_of_ratio = [[] for i in range(ClusterNum)]
        barcode_list_strain = [[] for i in range(ClusterNum + 1)]

        for i in range(snp_pd.shape[1]):
            column_name_temp = snp_pd.columns[i]
            if column_name_temp.startswith('strain') or column_name_temp.startswith('hierarchy_cluster'):
                continue

            if 'list' not in column_name_temp:
                y = column_name_temp.split('/')[-1]
            else:
                continue

            # if not column_name_temp in snp_pd.columns:
            # continue
            temp_df_check = snp_pd.loc[snp_pd[column_name_temp] != 0]
            for j in range(ClusterNum):
                temp_df_check_strain = temp_df_check.loc[temp_df_check['strain_' + str(j) + '_genotype'] != 0]
                total_count_same = sum(temp_df_check_strain[column_name_temp] == \
                                       temp_df_check_strain['strain_' + str(j) + '_genotype'])
                list_of_ratio[j].append(total_count_same / float(len(temp_df_check_strain[column_name_temp])))

            list_of_ratio_check = [list_of_ratio[j][-1] for j in range(ClusterNum)]
            print(list_of_ratio_check)

            if max(list_of_ratio_check) < 0.95:
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            elif sum(np.array(list_of_ratio_check) > 0.95) > 1:
                print('Not sure,more then two strains match', column_name_temp)
                # barcode_list_strain[-1].append(column_name_temp)
                barcode_list_strain[-1].append(y)
            else:
                # barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(column_name_temp)
                barcode_list_strain[list_of_ratio_check.index(max(list_of_ratio_check))].append(y)

        with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'w') as output:
            output.write('Cluster\tCell\n')

        for cluster in range(len(barcode_list_strain)):
            # print('look here is barcode_list_strain:', i)
            # print(barcode_list_strain[i])
            for cell in barcode_list_strain[cluster]:
                with open(os.path.join(PrepareDir, 'StrainCells.txt'), 'a') as output:
                    output.write(str(cluster) + '\t' + str(cell) + '\n')

    @timeit
    def AllBinSplit(self):
        # 需要输入想分菌株的ClusterxxxResult路径以及该箱判断分为几个簇。
        # BinClusterAnno.txt (SpeciesName必须与箱下的fasta文件名一致)
        # BinPrepareDir    ClusterNum    SpeciesName
        # /data_alluser/../SGB4933Result 4(分簇个数不能超过10)   SGB4933

        # 读取BinClusterAnno文件，调用SingleBinSplit()逐箱分菌株


        BinClusterAnnoDir = self.BinClusterAnno
        BinClusterAnno = pd.read_csv(BinClusterAnnoDir, sep='\t', header=0)
        for index, row in BinClusterAnno.iterrows():
            BinPrepareDir = row['BinPrepareDir']
            ClusterNum = int(row['ClusterNum'])
            SpeciesName = row['SpeciesName']

            #SingleBinSplit(BinPrepareDir, ClusterNum, SpeciesName)
            self.SingleBinSplit(BinPrepareDir,ClusterNum, SpeciesName)
