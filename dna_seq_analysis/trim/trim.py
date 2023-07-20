from pathlib import Path
import multiprocessing
from tqdm import tqdm
import unittest
import sys

from dna_seq_analysis.tools.common import *


TRIM_DICT = {'pe': {'LEADING':'3',
                    'TRAILING':'3',
                    'SLIDINGWINDOW':'4:15',
                    'MINLEN':'36'}
                    ,
             'se': {'LEADING':'3',
                    'TRAILING':'3',
                    'SLIDINGWINDOW':'4:15',
                    'MINLEN':'36'}
            }


# Trim reads
class Trim_reads():
    """
    trim reads
    """
    def __init__(self,fastq,outdir,wildcards,trimmer,mod,thread):
        self.fastq = fastq
        self.mod = mod
        self.outdir = Path(outdir)
        self.sample,self.unit = wildcards
        self.trimmer = trimmer
        self.threads = thread

        # log
        self.log_dir = self.outdir/"logs/trim_reads"
        self.log_dir.mkdir(parents=True,exist_ok=True)
    
    @staticmethod
    def distribute_threads(input_files, output_files, available_threads):
        """
        Distribute available threads between trimmomatic itself and any potential pigz instances
        """
        gzipped_input_files = sum(1 for file in input_files if file.endswith(".gz"))
        gzipped_output_files = sum(1 for file in output_files if file.endswith(".gz"))
        potential_threads_per_process = available_threads // (
            1 + gzipped_input_files + gzipped_output_files
        )
        if potential_threads_per_process > 0:
            # decompressing pigz creates at most 4 threads
            pigz_input_threads = (
                min(4, potential_threads_per_process) if gzipped_input_files != 0 else 0
            )
            pigz_output_threads = (
                (available_threads - pigz_input_threads * gzipped_input_files)
                // (1 + gzipped_output_files)
                if gzipped_output_files != 0
                else 0
            )
            trimmomatic_threads = (
                available_threads
                - pigz_input_threads * gzipped_input_files
                - pigz_output_threads * gzipped_output_files
            )
        else:
            # not enough threads for pigz
            pigz_input_threads = 0
            pigz_output_threads = 0
            trimmomatic_threads = available_threads
        return trimmomatic_threads, pigz_input_threads, pigz_output_threads

    @staticmethod
    def compose_input_gz(filename, threads):
        if filename.endswith(".gz") and threads > 0:
            return f'''<(pigz -p {threads} --decompress --stdout {filename})'''
        return filename

    @staticmethod
    def compose_output_gz(filename, threads, compression_level):
        if filename.endswith(".gz") and threads > 0:
            return f'''>(pigz -p {threads} {compression_level} > {filename})'''
        return filename


    def trim_read(self):
        _logs = self.log_dir/f'{self.sample}-{self.unit}.log'
        compression_level = "-5"
        # Distribute threads
        input_files = [self.fastq['r1'], self.fastq['r2']]
        trimmed_dir = self.outdir/"01.trimmed"
        trimmed_dir.mkdir(parents=True,exist_ok=True)
        output_files = [
            f'{trimmed_dir}/{self.sample}-{self.unit}.1.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.1.unpaired.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.2.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.2.unpaired.fastq.gz',
        ]

        trimmomatic_threads, input_threads, output_threads = self.distribute_threads(
            input_files, output_files, self.threads
        )
        # Abandonment
        input_r1, input_r2 = [
            self.compose_input_gz(filename, input_threads) for filename in input_files
        ]

        input_r1, input_r2 = self.fastq['r1'], self.fastq['r2']
        # Abandonment
        output_r1, output_r1_unp, output_r2, output_r2_unp = [
            self.compose_output_gz(filename, output_threads, compression_level)
            for filename in output_files
        ]
        
        output_r1, output_r1_unp, output_r2, output_r2_unp = output_files[0],output_files[1],output_files[2],output_files[3]
        
        cmd = (
                f"trimmomatic {self.mod} -threads {trimmomatic_threads} "
                f"-trimlog {trimmed_dir}/{self.sample}-{self.unit}.trimlog.txt "
                f"{input_r1} {input_r2} "
                f"{output_r1} {output_r1_unp} "
                f"{output_r2} {output_r2_unp} "
                f"{self.trimmer} > {_logs}"
            )
        debug_subprocess_call(cmd)


@add_log
def run(params):
    """
    """
    fastq,outdir,wildcards,trimmer,mod,thread = params
    app = Trim_reads(fastq=fastq,outdir=outdir,wildcards=wildcards,trimmer=trimmer,mod=mod,thread=thread)
    app.trim_read()


@add_log
def trim(args):
    config_path = args.config_path
    thread = args.thread
    trim_param = args.trim_param
    config = parse_config(config_path)
    outdir = config['outdir']
    units_file = f'{config_path}/units.tsv'
    
    extra_param_dict = {param.split(':',1)[0]:param.split(':',1)[1] for param in trim_param.split()}
    
    units = get_units(units_file)
    
    param_list = []
    for wildcards in units.index:
        fastq = get_fastq(units,wildcards)
        # Applicable to additional test data
        if len(fastq['r1'].split(",")) > 1:
            tmp_fastq = Path(f'{outdir}/tmp_fastq')
            tmp_fastq.mkdir(parents=True,exist_ok=True)
            fastq_dict = {}
            for key,value in fastq.items():
                if len(value) > 1:
                    fastq_list = value.split(",")
                    soft = "zcat" if fastq_list[0][-2:] == 'gz' else 'cat'
                    fastqs = " ".join(fastq_list)
                    cmd = f"{soft} {fastqs} > {tmp_fastq}/{wildcards[0]}_{key.upper()}.fastq"
                    debug_subprocess_call(cmd)
                    fastq_dict.update({key:f'{tmp_fastq}/{wildcards[0]}_{key.upper()}.fastq'})
            mod = 'PE' if len(fastq) > 1 else 'SE'
            TRIM_DICT[mod.swapcase()].update(extra_param_dict)
            trimmer = " ".join([key+':'+value for key,value in TRIM_DICT[mod.swapcase()].items()])
            param_list.append((fastq_dict,outdir,wildcards,trimmer,mod,thread)) 
        else:
            mod = 'PE' if len(fastq) > 1 else 'SE'
            TRIM_DICT[mod.swapcase()].update(extra_param_dict)
            trimmer = " ".join([key+':'+value for key,value in TRIM_DICT[mod.swapcase()].items()])
            param_list.append((fastq,outdir,wildcards,trimmer,mod,thread))
    
    with multiprocessing.Pool(len(param_list)) as p:
        list(tqdm(p.imap(run,param_list),total=len(param_list),unit_scale = True,ncols = 70,file = sys.stdout,desc='Trim reads '))
    p.close()
    p.join()
    
    # rm tmp_fastq
    tmp_fastq = Path(f'{outdir}/tmp_fastq')
    if tmp_fastq.exists():
        for fastq in tmp_fastq.rglob("*"):
            pathlib.Path(fastq).unlink()
        tmp_fastq.rmdir()

    
def get_opts_trim(parser, sub_program=True):
    parser.add_argument('--trim_param',help="Additional parameters for the called software. Need to be enclosed in quotation marks.\
For example, `--{software}_param 'param1:value1 param2:value2'`,`--trim_param 'LEADING:3 TRAILING:3'`.",default='LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')
    parser.add_argument('--thread',help='Number of threads.', default=4,type=int)
    
    if sub_program:
        parser = s_common(parser)
    return parser

if __name__ == '__main__':
    unittest.main()