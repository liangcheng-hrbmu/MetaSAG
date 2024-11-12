import os
import pandas as pd
import random
import subprocess
import time

'''
脚本注释
1. 本脚本需要输入一个FastaDir,还有一个树的中间文件路劲TreeTemp,返回一个树文件
2. 输入箱的注释文件 Bin Phylum CellNum用来作为itol绘制树的输入

'''
PHY_COLOR = {"Verrucomicrobiota":"#660000","Thermoproteota":"#783f04","Thermoplasmatota":"#7f6000","Spirochaetota":"#274e13","Proteobacteria":"#0c343d","Planctomycetota":"#073763","Patescibacteria":"#20124d","Myxococcota":"#4c1130","Marinisomatota":"#cc0000","Halobacteriota":"#e69138","Gemmatimonadota":"#f1c232","Firmicutes_C":"#6aa84f","Firmicutes_B":"#45818e","Firmicutes_A":"#3d85c6","Firmicutes":"#674ea7","Desulfobacterota_I":"#a64d79","Desulfobacterota":"#ea9999","Deinococcota":"#f9cb9c","Cyanobacteria":"#ffe599","Chloroflexota":"#b6d7a8","Campylobacterota":"#a2c4c9","Bacteroidota":"#9fc5e8","Actinobacteriota":"#b4a7d6","Acidobacteriota":"#d5a6bd","NA":"#000000","ERROR":"#ff0000"}



def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    return wrapper





def getRandomColor():
    red =  random.randint(0,255)
    green = random.randint(0,255)
    blue = random.randint(0,255)
    color_code = "#{:02x}{:02x}{:02x}".format(red,green,blue)
    return color_code

@timeit
def BuildTree(FastaDir,TreeTemp,env=None):
    #要求TreeTemp必须是绝对路径
    FastaAdjust = os.path.join(TreeTemp,'FastaAdjust/')
    DB = os.path.join(TreeTemp,'DB/')
    Anno = os.path.join(TreeTemp,'Anno/')

    # 创建文件夹
    DirCreat = [FastaAdjust, DB, Anno]
    for dir in DirCreat:
        os.makedirs(dir, exist_ok=True)

    #统计fasta文件
    files = os.listdir(FastaDir)
    fastafiles = []  #
    for filename in files:
        if filename.endswith('.fasta'):
            filename = filename[:-6]
            fastafiles.append(filename)

    # 对fasta构建数据库文件
    for fastafile in fastafiles:
        inputfasta = os.path.join(FastaDir,fastafile) + '.fasta'
        outputdb = os.path.join(DB,fastafile) + '.db'
        command = 'anvi-gen-contigs-database -f ' + inputfasta + ' -o ' + outputdb
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command,shell=True)

    # 构建不了的基因组文件
    for fastafile in fastafiles:
        inputfasta = os.path.join(FastaDir, fastafile) + '.fasta'
        outputdb = os.path.join(DB, fastafile) + '.db'
        adjustfasta = os.path.join(FastaAdjust, fastafile) + '.fasta'
        adjustfastaname = os.path.join(FastaAdjust, fastafile) + '_name_conversion.txt'
        if not os.path.exists(outputdb):
            command = 'anvi-script-reformat-fasta ' + inputfasta + ' -o ' + adjustfasta + ' --simplify --report ' + adjustfastaname
            if env is not None:
                subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)


            command = 'anvi-gen-contigs-database -f ' + adjustfasta +' -o ' + outputdb
            if env is not None:
                subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
            else:
                subprocess.call(command, shell=True)

    # 构建external.genomes.txt
    files = os.listdir(DB)
    DBfiles = []
    for filename in files:
        if filename.endswith('.db'):
            filename = filename[:-3]
            DBfiles.append(filename)

    externalfile = os.path.join(TreeTemp,'external_genomes.txt')
    concatenatedfile = os.path.join(TreeTemp,'concatenated-proteins_Bacteria_71_ribosomal6.fa')
    phylogenomicfile = os.path.join(TreeTemp, 'phylogenomic-tree_Bacteria_71_ribosomal6.txt')

    with open(externalfile, 'w') as input1:
        input1.write('name' + '\t' + 'contigs_db_path' + '\n')

    for db in DBfiles:
        with open(externalfile, 'a') as input1:
            input1.write(db + '\t' + os.path.join(DB,db) + '.db' + '\n')

    # 为每个db构建hmm模型
    for db in DBfiles:
        command = 'anvi-run-hmms -c ' + os.path.join(DB,db) + '.db'
        if env is not None:
            subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
        else:
            subprocess.call(command, shell=True)

    command = 'anvi-get-sequences-for-hmm-hits --external-genomes ' + externalfile + ' -o ' + concatenatedfile + ' --hmm-source Bacteria_71 --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 --return-best-hit --get-aa-sequences --concatenate'
    if env is not None:
        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)

    command = 'anvi-gen-phylogenomic-tree -f ' + concatenatedfile + ' -o ' + phylogenomicfile
    if env is not None:
        subprocess.run(['conda', 'run', '-n', env, 'bash', '-c', command])
    else:
        subprocess.call(command, shell=True)


#BinAnno文件格式
# Bin Phylum CellNum Color 第一种，为每个Phylum设置了Color(这种比较简单)
# Bin Phylum CellNum 第二种，没有为Phylum设置Color(自己为不常见Phylum随机设置颜色，常见Phylum的Color直接在PHY_COLOR中选择)

def itolPlot(BinAnno,Anno):
    #创建目录
    if not os.path.exists(Anno):
        os.makedirs(Anno,exist_ok=True)

    BinAnno = pd.read_csv(BinAnno,sep='\t',header=0)
    #创建两个itol_Anno文件：Label_Anno,bar_anno
    Label_Anno = os.path.join(Anno,'Label_Anno.txt')
    Bar_Anno = os.path.join(Anno,'Bar_Anno.txt')

    with open(Label_Anno, 'w') as input1:
        input1.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')

    with open(Bar_Anno, 'w') as input1:
        input1.write(
            'DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,simple bar testing\nDATASET_SCALE,100-1st line at 100-#5b5b5b-5-1-1,500-2nd line at 200-#5b5b5b-5-1-1,1000-3rd line at 1000-#5b5b5b-5-1-1,2000-4st line at 2000-#5b5b5b-5-1-1\nCOLOR,#999999\nWIDTH,1000\nMARGIN,0\nHEIGHT_FACTOR,1\nBAR_SHIFT,0\nBAR_ZERO,0\nDATA\n')

    #填写Label_Anno.txt
    if 'Color' in BinAnno.columns:
        for index, row in BinAnno.iterrows():
            phylum = row['Phylum']
            bin = row['Bin']
            color = row['Color']
            with open(Label_Anno,'a') as output:
                output.write(str(bin) + '\trange\t' + str(color) + '\t' + str(phylum) + '\n')
    else:
        for index, row in BinAnno.iterrows():

            phylum = row['Phylum']
            bin = row['Bin']

            if phylum in PHY_COLOR.keys():
                color = PHY_COLOR[phylum]
            else:
                #给非常见Phylum分配一个color并放到PHY_COLOR关键字中
                while True:
                    color = getRandomColor()
                    if color not in PHY_COLOR.values():
                        PHY_COLOR[phylum] = color
                        break
                color = PHY_COLOR[phylum]

            with open(Label_Anno, 'a') as output:
                output.write(str(bin) + '\trange\t' + str(color) + '\t' + str(phylum) + '\n')


    #填写Bar_Anno.txt
    for index, row in BinAnno.iterrows():
        with open(Bar_Anno,'a') as output:
            output.write(str(row['Bin'])+','+str(row['CellNum'])+'\n')
