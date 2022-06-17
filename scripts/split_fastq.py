import pysam
import pandas as pd
from xopen import xopen
import re
import os
from collections import defaultdict,OrderedDict
from tenacity import *
import pathlib
import subprocess
import argparse
import logging
import time
from datetime import timedelta
from functools import wraps



I = 0


def add_log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10)  # debug

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper


@add_log
def split_adapter(fa):
    """
    creat adapter_dict
    {"adapter1":xxxxxxx}
    """
    adapter_dict ={}
    with pysam.FastxFile(fa) as fh:
        for record in fh:
            header,seq = record.name, record.sequence
            adapter_dict[header] = seq
    return adapter_dict


@add_log
def parse_mapfile(mapfile):
    """
    
    """
    args = pd.read_table(mapfile).set_index("sample_name", drop=False)
    args_dict = args.to_dict('index')
    return args_dict


def hm_dis(s1, s2):
    dis = sum(e1 != e2 for e1, e2 in zip(s1, s2))
    return dis


def check_return_info(return_info):
    """
    check fun return
    """
    if return_info > 0:
        return True
    else:
        return False
    

class Split_fastq():
    """
    Split the original sequence and generate a Config file
    """
    def __init__(self,mapfile):
        self.mapfile = mapfile
        self.args_dict = parse_mapfile(self.mapfile)
    

    @add_log
    def run(self):
        """
        split adapter
        """
        args_dict = self.args_dict
        for sample in args_dict:
            sample_name = args_dict[sample]['sample_name']
            fa = args_dict[sample]['fa']
            fq1 = args_dict[sample]['fq1']
            fq2 = args_dict[sample]['fq2']
            fq_outdir = args_dict[sample]['fq_outdir']
            config_outdir = args_dict[sample]['config_outdir']
            mod = args_dict[sample]['mod']
        
            pathlib.Path(fq_outdir).mkdir(parents=True,exist_ok=True)
        
            adapter_dict = split_adapter(fa)
            adapter_dict_sorted = OrderedDict(sorted(adapter_dict.items(),key=lambda l: int(re.findall('\d+', l[0])[0])))
            
            index_len = len(adapter_dict_sorted)

            units_dict = defaultdict(dict)
            sample_dict = defaultdict(dict)
            for adapter in adapter_dict:
                sample_dict[f"{sample_name}_{adapter}"] = {"sample":f"{sample_name}_{adapter}"}
                units_dict[f"{sample_name}_{adapter}"] = {"sample":f"{sample_name}_{adapter}","unit":1,"platform":'ILLUMINA','fq1':f"{fq_outdir}/{sample_name}_{adapter}_R1.fastq",'fq2':f"{fq_outdir}/{sample_name}_{adapter}_R2.fastq"}
                
            sample_config = pd.DataFrame.from_dict(sample_dict).T
            sample_config.to_csv(f"{config_outdir}/{sample_name}_samples.tsv",sep="\t",index=None)
            units_config = pd.DataFrame.from_dict(units_dict).T
            units_config.to_csv(f"{config_outdir}/{sample_name}_units.tsv",sep="\t",index=None)
            
            
            self.adapter_list(sample_name,fq1,fq2,fq_outdir,mod,adapter_dict_sorted,index_len)
            

    @retry(retry=retry_if_result(check_return_info))
    def adapter_list(self,sample_name,fq1,fq2,fq_outdir,mod,adapter_dict_sorted,index_len):
        global I
        if mod == 'fq1':
            fq=fq1
        if mod == 'fq2':
            fq=fq2

        frist_adapter = list(adapter_dict_sorted.items())[0][0]
        name_list = xopen(f"{fq_outdir}/{sample_name}_{frist_adapter}.lst",'w')

        fh =  pysam.FastxFile(fq)
        for record in fh:
            header, seq, = record.name, record.sequence
            if hm_dis(seq[:9],adapter_dict_sorted[frist_adapter]) > 1:
                continue
            name_list.write(f'{header}\n')
        name_list.close()
        fh.close()

        out_fq_R1 = f'{fq_outdir}/{sample_name}_{frist_adapter}_R1.fastq'
        out_fq_R2 = f'{fq_outdir}/{sample_name}_{frist_adapter}_R2.fastq'

        cmd_line = (
                    f'seqtk subseq {fq1} {fq_outdir}/{sample_name}_{frist_adapter}.lst > {out_fq_R1};'
                    f'seqtk subseq {fq2} {fq_outdir}/{sample_name}_{frist_adapter}.lst > {out_fq_R2}'
                    )
       
        subprocess.check_call(cmd_line,shell=True)
        print(f'the {I+1} finished')
        
        if I < index_len:
            adapter_dict_sorted.popitem(last=False)
            I += 1
            i = index_len - I 
            return i
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile',help='mapfile',required=True)
    args = parser.parse_args()
    mapfile = args.mapfile
    runner =  Split_fastq(mapfile)
    runner.run()
    