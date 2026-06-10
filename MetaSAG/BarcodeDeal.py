import time
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns
import os
import shutil
import sys
import warnings
import subprocess
from itertools import islice
#Time Decorator
def timeit(func):
    def wrapper(*args,**kwargs):
        start_time = time.time()
        result = func(*args,**kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} took {elapsed_time:.4f} seconds to execute.")
        return result
    wrapper.__doc__=func.__doc__
    return wrapper




def _first_path(path_or_list):
    if isinstance(path_or_list, list):
        return path_or_list[0]
    return path_or_list



def _ensure_parent_dir(file_path):
    parent_dir = os.path.dirname(file_path)
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)



def _ensure_parent_dirs(file_paths):
    if isinstance(file_paths, list):
        for file_path in file_paths:
            _ensure_parent_dir(file_path)
    else:
        _ensure_parent_dir(file_paths)



def _remove_fastq(fastq):
    if isinstance(fastq, list):
        for file_path in fastq:
            os.remove(file_path)
    else:
        os.remove(fastq)



def _infer_reads_end(fastq):
    if isinstance(fastq, str):
        return 'Single'
    if isinstance(fastq, list):
        if len(fastq) == 1:
            return 'Single'
        if len(fastq) == 2:
            return 'Pair'
    raise ValueError("FASTQ input must be a string, a one-file list, or a two-file list.")



def _pair_paths(fastq):
    if not isinstance(fastq, list) or len(fastq) != 2:
        raise ValueError("Paired-end input must be a list of two FASTQ files.")

    if fastq[0].endswith('_R2.fastq') and fastq[1].endswith('_R1.fastq'):
        return fastq[1], fastq[0]
    return fastq[0], fastq[1]



def _infer_trim_prefix(InputFastq):
    ReadsEnd = _infer_reads_end(InputFastq)
    if ReadsEnd == 'Pair':
        InputFastq1, _ = _pair_paths(InputFastq)
        file_name = os.path.basename(InputFastq1)
        if file_name.endswith('_R1.fastq'):
            return file_name[:-len('_R1.fastq')]
    else:
        file_name = os.path.basename(_first_path(InputFastq))
        if file_name.endswith('.fastq'):
            return file_name[:-len('.fastq')]

    return os.path.splitext(file_name)[0]



def _cell_from_read_name(read_name):
    return read_name.split(':', 1)[0].replace('@', '')



def _remove_tag_from_read_name(read_name):
    if read_name.startswith('@') and ':' in read_name:
        return '@' + read_name.split(':', 1)[1]
    return read_name



def _parse_tag_result(tag_result):
    if not isinstance(tag_result, tuple) or len(tag_result) != 2:
        raise ValueError("TagFunc must return a tuple: (tag, status). Example: ('Cell500000', 'Valid') or ('BC', 'Error:BC').")

    ReadsTag, ReadsStatus = tag_result
    return str(ReadsTag), str(ReadsStatus)



def _format_reads_count(reads_count):
    reads_count = pd.DataFrame(
        [(status, tag, count) for (status, tag), count in reads_count.items()],
        columns=['Status', 'Tag', 'TagCount']
    )
    reads_count['StatusOrder'] = reads_count['Status'].apply(lambda status: 1 if status == 'Valid' else 0)
    reads_count = reads_count.sort_values(by=['StatusOrder', 'Status', 'Tag'])
    reads_count = reads_count.drop(columns=['StatusOrder'])
    reads_count = reads_count.set_index(['Status', 'Tag'])
    return reads_count



def _save_reads_count(reads_count, reads_count_file):
    if reads_count_file != 'No':
        _ensure_parent_dir(reads_count_file)
        reads_count.to_csv(reads_count_file)



@timeit
def BarcodeTag(InputFastq,CellFastq,ErrorFastq,TagFunc,RemoveInput=False,ReadsCountFile='No',ChunkSize=10000):

    ReadsCount = {} # Count reads by barcode type
    ReadsEnd = _infer_reads_end(InputFastq)

    if ReadsEnd=='Single':
        InputFastq = _first_path(InputFastq)
        CellFastq = _first_path(CellFastq)
        ErrorFastq = _first_path(ErrorFastq)
        _ensure_parent_dirs([CellFastq, ErrorFastq])

        with open(InputFastq,'r') as fin, open(CellFastq,'w') as cell_out, open(ErrorFastq,'w') as error_out:
            while True:
                lines=list(islice(fin,4))
                if not lines:
                    break
                if len(lines) % 4 != 0:
                    raise ValueError("Read sequence line count ≠ 4.")

                ReadsName = lines[0]
                ReadsSeq = lines[1]
                ReadsQC = lines[3]
                ReadsTag, ReadsStatus = _parse_tag_result(TagFunc(ReadsSeq))

                ReadsName = "@" + ReadsTag + ":" + ReadsName[1:] # Add Tag to ReadsName
                ReadsCountKey = (ReadsStatus, ReadsTag)

                if ReadsCountKey in ReadsCount.keys():
                    ReadsCount[ReadsCountKey] = ReadsCount[ReadsCountKey] + 1
                else:
                    ReadsCount[ReadsCountKey] = 1

                if ReadsStatus == "Valid":
                    cell_out.write(ReadsName + ReadsSeq + "+\n" + ReadsQC)
                else:
                    error_out.write(ReadsName + ReadsSeq + "+\n" + ReadsQC)

        if RemoveInput:
            _remove_fastq(InputFastq)

        ReadsCount = _format_reads_count(ReadsCount)
        _save_reads_count(ReadsCount, ReadsCountFile)
        return ReadsCount
    elif ReadsEnd=='Pair':

        InputFastq1, InputFastq2 = _pair_paths(InputFastq)

        CellFastq1 = CellFastq[0]
        CellFastq2 = CellFastq[1]

        ErrorFastq1 = ErrorFastq[0]
        ErrorFastq2 = ErrorFastq[1]
        _ensure_parent_dirs([CellFastq1, CellFastq2, ErrorFastq1, ErrorFastq2])

        with open(InputFastq1, 'r') as fin1, open(InputFastq2, 'r') as fin2, open(CellFastq1, 'w') as cell_out1, open(CellFastq2, 'w') as cell_out2, open(ErrorFastq1, 'w') as error_out1, open(ErrorFastq2, 'w') as error_out2:
            while True:
                lines1 = list(islice(fin1, 4))
                lines2 = list(islice(fin2, 4))
                if not lines1 and not lines2:
                    break
                if len(lines1) != 4 or len(lines2) != 4:
                    raise ValueError("Read sequence line count ≠ 4.")

                ReadsName1 = lines1[0]
                ReadsSeq1 = lines1[1]
                ReadsQC1 = lines1[3]
                ReadsTag, ReadsStatus = _parse_tag_result(TagFunc(ReadsSeq1))
                ReadsName1 = '@' + ReadsTag + ":" + ReadsName1[1:]  # Add Tag to ReadsName

                ReadsName2 = lines2[0]
                ReadsSeq2 = lines2[1]
                ReadsQC2 = lines2[3]
                ReadsName2 = '@' + ReadsTag + ":" + ReadsName2[1:]  # Add Tag to ReadsName
                ReadsCountKey = (ReadsStatus, ReadsTag)

                if ReadsCountKey in ReadsCount.keys():
                    ReadsCount[ReadsCountKey] = ReadsCount[ReadsCountKey] + 1
                else:
                    ReadsCount[ReadsCountKey] = 1

                if ReadsStatus == "Valid":
                    cell_out1.write(ReadsName1 + ReadsSeq1 + "+\n" + ReadsQC1)
                    cell_out2.write(ReadsName2 + ReadsSeq2 + "+\n" + ReadsQC2)
                else:
                    error_out1.write(ReadsName1 + ReadsSeq1 + "+\n" + ReadsQC1)
                    error_out2.write(ReadsName2 + ReadsSeq2 + "+\n" + ReadsQC2)

        if RemoveInput:
            _remove_fastq([InputFastq1, InputFastq2])

        ReadsCount = _format_reads_count(ReadsCount)
        _save_reads_count(ReadsCount, ReadsCountFile)
        return ReadsCount


    else:
        raise ValueError("Cannot infer whether input FASTQ is single-end or paired-end.")



PYTHONDIR = os.path.dirname(os.path.abspath(__file__))
toolDir=os.path.join(PYTHONDIR,"trimmomatic")
@timeit
def trim(InputFastq,TrimDir,TrimPrefix='No',ILLUMINACLIP='TruSeq3-PE.fa:2:30:10:3:TRUE',LEADING=25,TRAILING=3,SLIDINGWINDOW='4:20',MINLEN=30,threads=12,RemoveInput=False):

    # Create directory for trimmed results.
    if not os.path.exists(TrimDir):
        os.makedirs(TrimDir, exist_ok=True)

    ReadsEnd = _infer_reads_end(InputFastq)
    if TrimPrefix == 'No':
        TrimPrefix = _infer_trim_prefix(InputFastq)

    if ReadsEnd == 'Single':

        InputFastq = _first_path(InputFastq)

        TrimFastq = os.path.join(TrimDir, TrimPrefix) + '.fastq'
        TrimLog = os.path.join(TrimDir, TrimPrefix) + '.log'
        TrimError = os.path.join(TrimDir, TrimPrefix) + '.trim.stderr.txt'
        command = f"(java -jar {toolDir}/trimmomatic-0.36.jar SE -threads {threads} -phred33 \
                         -trimlog {TrimLog} \
                          {InputFastq} {TrimFastq} \
                           ILLUMINACLIP:{toolDir}/{ILLUMINACLIP} \
                            LEADING:{LEADING} \
                             TRAILING:{TRAILING} \
                              SLIDINGWINDOW:{SLIDINGWINDOW} \
                               MINLEN:{MINLEN}) \
                                2> {TrimError}"
        return_code = subprocess.call(command, shell=True)
        if return_code != 0:
            raise RuntimeError(f"Trimmomatic failed with exit code {return_code}. See {TrimError}.")


    elif ReadsEnd == 'Pair':
        InputFastq1, InputFastq2 = _pair_paths(InputFastq)

        TrimFastq1 = os.path.join(TrimDir, TrimPrefix) + '_R1.fastq'
        TrimFastq2 = os.path.join(TrimDir, TrimPrefix) + '_R2.fastq'
        TrimLog = os.path.join(TrimDir, TrimPrefix) + '.log'
        TrimUnpaired1 = os.path.join(TrimDir, TrimPrefix) + '_unpaired_R1.fastq'
        TrimUnpaired2 = os.path.join(TrimDir, TrimPrefix) + '_unpaired_R2.fastq'
        TrimError = os.path.join(TrimDir, TrimPrefix) + '.trim.stderr.txt'
        command = f"(java -jar {toolDir}/trimmomatic-0.36.jar PE -threads {threads} -phred33 \
             -trimlog {TrimLog} \
              {InputFastq1} {InputFastq2} {TrimFastq1} {TrimUnpaired1}  {TrimFastq2} {TrimUnpaired2} \
               ILLUMINACLIP:{toolDir}/{ILLUMINACLIP} \
                LEADING:{LEADING} \
                 TRAILING:{TRAILING} \
                  SLIDINGWINDOW:{SLIDINGWINDOW} \
                   MINLEN:{MINLEN}) \
                    2> {TrimError}"
        return_code = subprocess.call(command, shell=True)
        if return_code != 0:
            raise RuntimeError(f"Trimmomatic failed with exit code {return_code}. See {TrimError}.")

    if RemoveInput:
        _remove_fastq(InputFastq)

@timeit
def SplitTaggedFastqByCell(CellFastq,CellBarnDir,CellList:list = [],RemoveInput=False):
    if not os.path.exists(CellBarnDir):
        os.makedirs(CellBarnDir, exist_ok=True)

    ReadsEnd = _infer_reads_end(CellFastq)

    if CellList:
        if ReadsEnd == 'Single':
            # 创建所有细胞的fastq空文件
            for Cell in CellList:
                with open(os.path.join(CellBarnDir,Cell) + '.fastq','w') as f:
                    pass

            CellFastq = _first_path(CellFastq)

            with open(CellFastq, 'r') as f:
                while True:
                    lines = list(islice(f, 4))
                    if not lines:
                        break
                    if len(lines) != 4:
                        raise ValueError("Read sequence line count ≠ 4.")

                    bc = _cell_from_read_name(lines[0])

                    if bc in CellList:
                        lines[0] = _remove_tag_from_read_name(lines[0])
                        fastq = "".join(lines)
                        with open(os.path.join(CellBarnDir, str(bc)) + '.fastq', 'a') as fout:
                            fout.write(fastq)

        elif ReadsEnd == 'Pair':
            # 创建所有细胞的fastq空文件
            for Cell in CellList:
                with open(os.path.join(CellBarnDir, Cell) + '_R1.fastq', 'w') as f1, open(os.path.join(CellBarnDir, Cell) + '_R2.fastq','w') as f2:
                    pass

            CellFastq1, CellFastq2 = _pair_paths(CellFastq)

            with open(CellFastq1, 'r') as f1, open(CellFastq2, 'r') as f2:
                while True:
                    lines1 = list(islice(f1, 4))
                    lines2 = list(islice(f2, 4))
                    if not lines1 and not lines2:
                        break
                    if len(lines1) != 4 or len(lines2) != 4:
                        raise ValueError("Read sequence line count ≠ 4.")

                    bc = _cell_from_read_name(lines1[0])

                    if bc in CellList:
                        lines1[0] = _remove_tag_from_read_name(lines1[0])
                        lines2[0] = _remove_tag_from_read_name(lines2[0])
                        fastq1 = "".join(lines1)
                        fastq2 = "".join(lines2)
                        with open(os.path.join(CellBarnDir, str(bc)) + '_R1.fastq', 'a') as fout1, open(os.path.join(CellBarnDir, str(bc)) + '_R2.fastq','a') as fout2:
                            fout1.write(fastq1)
                            fout2.write(fastq2)

    else:
        if ReadsEnd == 'Single':
            CellFastq = _first_path(CellFastq)

            with open(CellFastq, 'r') as f:
                while True:
                    lines = list(islice(f, 4))
                    if not lines:
                        break
                    if len(lines) != 4:
                        raise ValueError("Read sequence line count ≠ 4.")

                    bc = _cell_from_read_name(lines[0])

                    lines[0] = _remove_tag_from_read_name(lines[0])
                    fastq = "".join(lines)
                    with open(os.path.join(CellBarnDir, str(bc)) + '.fastq', 'a') as fout:
                        fout.write(fastq)


        elif ReadsEnd == 'Pair':
            CellFastq1, CellFastq2 = _pair_paths(CellFastq)

            with open(CellFastq1, 'r') as f1, open(CellFastq2, 'r') as f2:
                while True:
                    lines1 = list(islice(f1, 4))
                    lines2 = list(islice(f2, 4))
                    if not lines1 and not lines2:
                        break
                    if len(lines1) != 4 or len(lines2) != 4:
                        raise ValueError("Read sequence line count ≠ 4.")

                    bc = _cell_from_read_name(lines1[0])

                    lines1[0] = _remove_tag_from_read_name(lines1[0])
                    lines2[0] = _remove_tag_from_read_name(lines2[0])
                    fastq1 = "".join(lines1)
                    fastq2 = "".join(lines2)
                    with open(os.path.join(CellBarnDir, str(bc)) + '_R1.fastq', 'a') as fout1, open(os.path.join(CellBarnDir, str(bc)) + '_R2.fastq','a') as fout2:
                        fout1.write(fastq1)
                        fout2.write(fastq2)

    if RemoveInput:
        _remove_fastq(CellFastq)


#CellDirt 是一个两列数据框，第一列为Cell名称，第二列为Cell序列文件的路径名如 /DATA/DMY/XXX/XXX/Cell123 (注意没有.fastq或_R1.fastq,_R2.fastq)
@timeit
def SplitTaggedFastqByCellPath(CellFastq,CellDirt,RemoveInput=False):
    CellDirt.columns = ['Cell','Dirt']
    CellDirt = CellDirt.drop_duplicates(subset=['Cell'],keep='first')
    CellDirt = dict(zip(CellDirt['Cell'],CellDirt['Dirt']))
    CellList = CellDirt.keys()
    ReadsEnd = _infer_reads_end(CellFastq)


    if ReadsEnd == 'Single':
        CellFastq = _first_path(CellFastq)

        with open(CellFastq, 'r') as f:
            while True:
                lines = list(islice(f, 4))
                if not lines:
                    break
                if len(lines) != 4:
                    raise ValueError("Read sequence line count ≠ 4.")

                bc = _cell_from_read_name(lines[0])

                if bc in CellList:
                    lines[0] = _remove_tag_from_read_name(lines[0])
                    fastq = "".join(lines)
                    file = CellDirt[bc]
                    file_dir = os.path.dirname(file)
                    if file_dir and not os.path.exists(file_dir):
                        os.makedirs(file_dir, exist_ok=True)
                    with open(file + '.fastq', 'a') as fout:
                        fout.write(fastq)


    elif ReadsEnd == 'Pair':

        CellFastq1 = CellFastq[0]
        CellFastq2 = CellFastq[1]

        with open(CellFastq1, 'r') as f1, open(CellFastq2, 'r') as f2:
            while True:
                lines1 = list(islice(f1, 4))
                lines2 = list(islice(f2, 4))
                if not lines1 and not lines2:
                    break
                if len(lines1) != 4 or len(lines2) != 4:
                    raise ValueError("Read sequence line count ≠ 4.")

                bc = _cell_from_read_name(lines1[0])

                if bc in CellList:
                    lines1[0] = _remove_tag_from_read_name(lines1[0])
                    lines2[0] = _remove_tag_from_read_name(lines2[0])
                    fastq1 = "".join(lines1)
                    fastq2 = "".join(lines2)
                    file = CellDirt[bc]
                    file_dir = os.path.dirname(file)
                    if file_dir and not os.path.exists(file_dir):
                        os.makedirs(file_dir, exist_ok=True)
                    with open(file + '_R1.fastq', 'a') as fout1, open(file + '_R2.fastq', 'a') as fout2:
                        fout1.write(fastq1)
                        fout2.write(fastq2)

    if RemoveInput:
        _remove_fastq(CellFastq)


class Barcode:
    def __init__(self, InputFastq, resultDir='No', SamplePrefix='No'):
        self.InputFastq = InputFastq
        self.resultDir = resultDir
        self.SamplePrefix = _infer_trim_prefix(InputFastq) if SamplePrefix == 'No' else SamplePrefix
        self.ReadsEnd = _infer_reads_end(InputFastq)

        self.TempDir = self._default_dir('temp')
        self.TrimDir = self._default_dir('Trim')
        self.CellBarnDir = self._default_dir('Cell')
        self.ReadsCountFile = self._default_file('ReadsCount.csv')

        self.CellFastq = self._default_fastq(self.TempDir, self.SamplePrefix)
        self.ErrorFastq = self._default_fastq(self.TempDir, f'{self.SamplePrefix}_error')
        self.TrimFastq = self._default_fastq(self.TrimDir, self.SamplePrefix)

    def BarcodeTag(self, CellFastq='No', ErrorFastq='No', TagFunc='No', RemoveInput=False, ReadsCountFile='No', ChunkSize=10000):
        if TagFunc == 'No':
            raise ValueError('TagFunc is required for Barcode.BarcodeTag().')

        CellFastq = self._resolve_path(CellFastq, self.CellFastq, 'CellFastq')
        ErrorFastq = self._resolve_path(ErrorFastq, self.ErrorFastq, 'ErrorFastq')
        ReadsCountFile = self._resolve_optional_path(ReadsCountFile, self.ReadsCountFile)

        reads_count = globals()['BarcodeTag'](
            InputFastq=self.InputFastq,
            CellFastq=CellFastq,
            ErrorFastq=ErrorFastq,
            TagFunc=TagFunc,
            RemoveInput=RemoveInput,
            ReadsCountFile=ReadsCountFile,
            ChunkSize=ChunkSize
        )

        self.CellFastq = CellFastq
        self.ErrorFastq = ErrorFastq
        return reads_count

    def trim(self, InputFastq='No', TrimDir='No', TrimPrefix='No', ILLUMINACLIP='TruSeq3-PE.fa:2:30:10:3:TRUE', LEADING=25, TRAILING=3, SLIDINGWINDOW='4:20', MINLEN=30, threads=12, RemoveInput=False):
        InputFastq = self._resolve_existing_fastq(InputFastq, [self.CellFastq, self.InputFastq], 'InputFastq')
        TrimDir = self._resolve_path(TrimDir, self.TrimDir, 'TrimDir')

        trim(
            InputFastq=InputFastq,
            TrimDir=TrimDir,
            TrimPrefix=TrimPrefix,
            ILLUMINACLIP=ILLUMINACLIP,
            LEADING=LEADING,
            TRAILING=TRAILING,
            SLIDINGWINDOW=SLIDINGWINDOW,
            MINLEN=MINLEN,
            threads=threads,
            RemoveInput=RemoveInput
        )

        self.TrimDir = TrimDir
        self.TrimFastq = self._default_fastq(TrimDir, _infer_trim_prefix(InputFastq) if TrimPrefix == 'No' else TrimPrefix)
        return self.TrimFastq

    def SplitTaggedFastqByCell(self, CellFastq='No', CellBarnDir='No', CellList:list = [], RemoveInput=False):
        CellFastq = self._resolve_existing_fastq(CellFastq, [self.TrimFastq, self.CellFastq], 'CellFastq')
        CellBarnDir = self._resolve_path(CellBarnDir, self.CellBarnDir, 'CellBarnDir')

        SplitTaggedFastqByCell(
            CellFastq=CellFastq,
            CellBarnDir=CellBarnDir,
            CellList=CellList,
            RemoveInput=RemoveInput
        )

        self.CellBarnDir = CellBarnDir

    def SplitTaggedFastqByCellPath(self, CellFastq='No', CellDirt=None, RemoveInput=False):
        if CellDirt is None:
            raise ValueError('CellDirt is required for Barcode.SplitTaggedFastqByCellPath().')

        CellFastq = self._resolve_existing_fastq(CellFastq, [self.TrimFastq, self.CellFastq], 'CellFastq')

        SplitTaggedFastqByCellPath(
            CellFastq=CellFastq,
            CellDirt=CellDirt,
            RemoveInput=RemoveInput
        )

    def _default_dir(self, dirname):
        if self.resultDir == 'No':
            return 'No'
        return os.path.join(self.resultDir, dirname)

    def _default_file(self, filename):
        if self.resultDir == 'No':
            return 'No'
        return os.path.join(self.resultDir, filename)

    def _default_fastq(self, directory, prefix):
        if directory == 'No':
            return 'No'
        if self.ReadsEnd == 'Pair':
            return [
                os.path.join(directory, f'{prefix}_R1.fastq'),
                os.path.join(directory, f'{prefix}_R2.fastq')
            ]
        return os.path.join(directory, f'{prefix}.fastq')

    def _resolve_path(self, user_path, default_path, name):
        if user_path != 'No':
            return user_path
        if default_path != 'No':
            return default_path
        raise ValueError(f'{name} is required when resultDir is not provided.')

    def _resolve_optional_path(self, user_path, default_path):
        if user_path != 'No':
            return user_path
        return default_path

    def _resolve_existing_fastq(self, user_path, candidate_paths, name):
        if user_path != 'No':
            return user_path

        for candidate_path in candidate_paths:
            if candidate_path != 'No' and self._fastq_exists(candidate_path):
                return candidate_path

        for candidate_path in candidate_paths:
            if candidate_path != 'No':
                return candidate_path

        raise ValueError(f'{name} is required when resultDir is not provided.')

    def _fastq_exists(self, fastq):
        if isinstance(fastq, list):
            return all(os.path.exists(file) for file in fastq)
        return os.path.exists(fastq)
