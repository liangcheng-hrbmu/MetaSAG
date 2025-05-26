import os
import pandas as pd
import subprocess
import time

'''

脚本注释
这个脚本的目的是：
1. [Diamond] 对于输入的fastq文件[fastq文件的介绍在后面会有不同类型的介绍],首先进行Diamond比对，得到fastq文件映射的Uniref片段注释。将Diamond结果文件整理得到Cell-Uniref Count矩阵；
2. [Seurat] 对于输入的Cell-Uniref Count矩阵，使用SeuratSGB.R脚本进行细胞的聚类，以及细胞的SGB - Cluster统计绘图（绘图部分暂未实现）
3. [HUMAnN] 整理出每一个Cluster细胞的Diamond注释结果，进行HUMAnN的通路注释，并使用PathClusterHeatmap.R绘制Cluster - Path的coverage_rate热图。

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



class HP():
    def __init__(self,FastqDir,ResultDir):
        self.FastqDir = FastqDir
        self.ResultDir = ResultDir
        if not os.path.exists(self.ResultDir):
            os.makedirs(self.ResultDir,exist_ok=True)

    PYTHONDIR = os.path.dirname(os.path.abspath(__file__))

    @timeit
    def Diamond(self, DiamondDB=None, Diamondenv=None):

        FastqDir = self.FastqDir
        DiamondDir = os.path.join(self.ResultDir,'DiamondDir')
        self.DiamondDir = DiamondDir

        # 本函数所输入的Fastq文件可以为一个，也可以是多个。

        # 创建目录
        if not os.path.exists(DiamondDir):
            os.makedirs(DiamondDir, exist_ok=True)

        # 遍历fastq文件进行Diamond得到xxx_uniref
        files = os.listdir(FastqDir)
        for filename in files:
            if filename.endswith('.fastq'):
                fastqfile = os.path.join(FastqDir, filename)
                DiamondFile = os.path.join(DiamondDir, filename[:-6] + '_uniref')
                if DiamondDB is not None:
                    command = 'diamond blastx --db ' + DiamondDB + ' -q ' + fastqfile + ' -o ' + DiamondFile + ' -k 1'
                else:
                    command = 'diamond blastx -q ' + fastqfile + ' -o ' + DiamondFile + ' -k 1'

                if Diamondenv is not None:
                    subprocess.run(['conda', 'run', '-n', Diamondenv, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)


    @timeit
    def Uniref2Matrix(self, MinUnirefNum=5, MinCellNum=5):
        DiamondDir = self.DiamondDir
        MatrixDir = os.path.join(self.ResultDir,'MatrixDir')
        self.MatrixDir = MatrixDir

        # 创建目录
        if not os.path.exists(MatrixDir):
            os.makedirs(MatrixDir)

        # 输入比对得到的所有Uniref结果，整理为count cell uniref file（cell所属的diamond_uniref文件名）
        AllSGBCellUnirefCount = pd.DataFrame(columns=['count', 'cell', 'uniref', 'file'])
        for unirefFile in os.listdir(DiamondDir):
            if unirefFile.endswith('_uniref'):
                unirefFileDir = os.path.join(DiamondDir, unirefFile)

                CellUnirefCount = {}
                with open(unirefFileDir) as input1:
                    while True:
                        line = input1.readline()
                        if len(line) == 0:
                            break

                        items = line.split('\t')
                        cell = items[0].split(':')[0]
                        uniref = items[1].split('|')[0]
                        ident = float(items[2])

                        key = str(cell) + '@' + str(uniref)
                        if ident > 90:
                            if key in CellUnirefCount.keys():
                                CellUnirefCount[key] += 1
                            else:
                                CellUnirefCount[key] = 1

                CellUnirefCount = pd.DataFrame(list(CellUnirefCount.items()))
                df = CellUnirefCount[0].str.split('@', expand=True)
                CellUnirefCount.columns = ['cell_uniref', 'count']
                CellUnirefCount["cell"] = df[0]
                CellUnirefCount['uniref'] = df[1]
                CellUnirefCount = CellUnirefCount.drop(columns='cell_uniref')
                CellUnirefCount['file'] = unirefFile
                CellUnirefCount.to_csv(os.path.join(MatrixDir, unirefFile[:-7] + '_CellUnirefCount.txt'), sep='\t',
                                       index=False)
                AllSGBCellUnirefCount = pd.concat([AllSGBCellUnirefCount, CellUnirefCount], ignore_index=True)

        AllSGBCellUnirefCount.to_csv(os.path.join(MatrixDir, 'AllSGBCellUnirefCount.txt'), sep='\t', index=False)
        self.AllSGBCellUnirefCount = AllSGBCellUnirefCount

        # 将所有细胞整理为cell-uniref count矩阵并输出文件
        # 去除在小于5个细胞中出现的Uniref90
        data = AllSGBCellUnirefCount
        Uniref10 = data['uniref'].value_counts()
        Uniref10Filter = [uniref for uniref in Uniref10.index if Uniref10[uniref] < MinCellNum]
        data = data[~data['uniref'].isin(Uniref10Filter)]
        # 去除存在少于5个uniref的细胞
        Cell10 = data['cell'].value_counts()
        Cell10Filter = [cell for cell in Cell10.index if Cell10[cell] < MinUnirefNum]
        data = data[~data['cell'].isin(Cell10Filter)]
        matrix = data.pivot(index='uniref', columns='cell', values='count').fillna(0)
        self.AllSGB_matrix = os.path.join(MatrixDir, 'AllSGB_matrix')
        matrix.to_csv(os.path.join(MatrixDir, 'AllSGB_matrix'), sep='\t')

    # 输入AllSGB_matrix进行细胞聚类，得到聚类结果
    @timeit
    def SeuratCluster(self,nfeatures=20000,npcs=30,FindNeighborsDim="1:10",FindClustersRes=0.5,UmapPDFWidth=6,UmapPDFHeight=4,DimPlotMethod='umap',FindAllMarkersPCT=0.25,FindAllMarkerslogFC=0.25,topN=10,HeatmapPDFWidth=12,HeatmapPDFHeight=8):

        SeuratResult = os.path.join(self.ResultDir,'SeuratResult')
        AllSGB_matrix = self.AllSGB_matrix
        # 创建目录
        if not os.path.exists(SeuratResult):
            os.makedirs(SeuratResult, exist_ok=True)
        SeuratSGB_R = os.path.join(self.PYTHONDIR, 'SeuratSGB.R')
        #command = 'Rscript ' + SeuratSGB_R + ' -i ' + AllSGB_matrix + ' -o ' + SeuratResult
        command =f'Rscript {SeuratSGB_R} -i {AllSGB_matrix} -o {SeuratResult} -nfeatures {nfeatures} -npcs {npcs} -FindNeighborsDim {FindNeighborsDim} -FindClustersRes {FindClustersRes} -UmapPDFWidth {UmapPDFWidth} -UmapPDFHeight {UmapPDFHeight} -DimPlotMethod {DimPlotMethod} -FindAllMarkersPCT {FindAllMarkersPCT} -FindAllMarkerslogFC {FindAllMarkerslogFC} -topN {topN} -HeatmapPDFWidth {HeatmapPDFWidth} -HeatmapPDFHeight {HeatmapPDFHeight}'
        subprocess.call(command, shell=True)

    @timeit
    def getPathCovCluster(self,HUMAnNFolderDirt):
        HUMAnNFolders = [folder for folder in os.listdir(HUMAnNFolderDirt) if folder.endswith('_humann')]

        PathCovCluster = pd.DataFrame(columns=['Path', 'Coverage', 'Cluster'])
        for folder in HUMAnNFolders:
            cluster = folder.replace('_humann', '')
            for file_name in os.listdir(os.path.join(HUMAnNFolderDirt, folder)):
                if file_name.endswith('_pathcoverage.tsv'):
                    file_path = os.path.join(HUMAnNFolderDirt, folder, file_name)
                    temp = pd.read_csv(file_path, sep='\t', header=None, comment='#')
                    temp.columns = ['Path', 'Coverage']
                    temp = temp[(~temp['Path'].str.contains('\|')) & (temp['Coverage'] > 0.1)]
                    # temp = temp[(~temp['Path'].str.contains('UNMAPPED')) | (~temp['Path'].str.contains('UNINTEGRATED'))]
                    temp = temp[~temp['Path'].str.contains('UNMAPPED|UNINTEGRATED')]
                    temp['Cluster'] = cluster

                    PathCovCluster = pd.concat([PathCovCluster, temp], ignore_index=True)
        PathCovCluster.to_csv(os.path.join(HUMAnNFolderDirt, 'PathCovCluster.txt'), sep='\t', index=False)

    # 分组细胞的uniref文件：
    # 输入CellAnno文件，请再输入细胞分组的列名（eg. “SGB”,"Cluster",“Sam”,“DiseaseState”,推荐使用“SGB”或“Cluster”）
    # 要求输入的文件：CellAnno
    @timeit
    def HUMAnNPath(self, CellAnno, Group, HUMAnNenv=None):

        DiamondDir = self.DiamondDir
        HUMAnNResult = os.path.join(self.ResultDir,'HUMAnNResult')

        UnirefForHUMAnN = os.path.join(HUMAnNResult, 'UnirefForHUMAnN')
        HUMAnNPathForGroup = os.path.join(HUMAnNResult, 'HUMAnNPathForGroup')
        HeatMapPlot = os.path.join(HUMAnNResult, 'HeatMapPlot')

        # 创建目录
        DirCreate = [HUMAnNResult, UnirefForHUMAnN, HUMAnNPathForGroup, HeatMapPlot]
        for dir in DirCreate:
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True)

        # 检索DiamondDir下所有的uniref文件，根据CellAnno文件以及group将uniref分到相应group的UnirefForHUMAnN/group_uniref文件下
        CellAnno = pd.read_csv(CellAnno, sep='\t', header=0)
        Cell_Group = dict(zip(CellAnno["Cell"], CellAnno[Group]))

        for file in os.listdir(DiamondDir):
            if file.endswith('_uniref'):
                UnirefFile = os.path.join(DiamondDir, file)
                with open(UnirefFile) as input1:
                    while True:
                        line = input1.readline()
                        if len(line) == 0:
                            break

                        cell = line.split(':')[0]
                        group = Cell_Group[cell]

                        groupUnirefFile = os.path.join(UnirefForHUMAnN, str(group) + '_uniref')
                        with open(groupUnirefFile, 'a') as output:
                            output.write(line)

        # 对UnirefForHUMAnN中每个uniref文件进行HUMAnN通路注释，结果文件放在HUMAnNPathForGroup文件夹下
        files = os.listdir(UnirefForHUMAnN)
        for filename in files:
            if filename.endswith('_uniref'):
                unirefFileDir = os.path.join(UnirefForHUMAnN, filename)
                humannFileDir = os.path.join(HUMAnNPathForGroup, filename.replace('_uniref', '')) + '_humann'
                command = 'humann --input ' + unirefFileDir + ' --output ' + humannFileDir
                if HUMAnNenv is not None:
                    subprocess.run(['conda', 'run', '-n', HUMAnNenv, 'bash', '-c', command])
                else:
                    subprocess.call(command, shell=True)

        # 整合HUMAnNPathForGroup文件夹下所有的group_humann文件
        self.getPathCovCluster(HUMAnNPathForGroup)
        # 绘制Heatmap热图,并将热图相关结果放在HeatMapPlot文件夹下
        PathCovCluster = os.path.join(HUMAnNPathForGroup, 'PathCovCluster.txt')
        PathClusterHeatmap_R = os.path.join(self.PYTHONDIR, 'PathClusterHeatmap.R')
        command = 'Rscript ' + PathClusterHeatmap_R + ' -i ' + PathCovCluster + ' -o ' + HeatMapPlot + '/'
        subprocess.call(command, shell=True)

