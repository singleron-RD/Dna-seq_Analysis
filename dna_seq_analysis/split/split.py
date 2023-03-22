import pandas as pd
import pysam
from pathlib import Path
import multiprocessing
from tqdm import tqdm
from collections import defaultdict
from xopen import xopen
import unittest

from dna_seq_analysis.tools.common import *


def split_adapter(fa):
    """
    creat adapter_dict
    {"adapter1":xxxxxxx}
    """
    adapter_dict = {}
    with pysam.FastxFile(fa) as fh:
        for record in fh:
            header,seq = record.name, record.sequence
            adapter_dict[header] = seq
    return adapter_dict


def parse_mapfile(mapfile):
    """
    
    """
    args = pd.read_table(mapfile).set_index("sample_name", drop=False)
    args_dict = args.to_dict('index')
    return args_dict


def hm_dis(s1, s2):
    dis = sum(e1 != e2 for e1, e2 in zip(s1, s2))
    return dis


class Split_fastq():
    """
    Split the original sequence and generate a Config file
    """
    def __init__(self,args, display_title=None):
        
        Step.__init__(self, args, display_title=display_title)
        
        self.config_path = args.config_path
        self.config = parse_config(self.config_path)
        self.mapfile = args.mapfile
        self.outdir = self.config['outdir']
        self.args_dict = parse_mapfile(self.mapfile)
        self.units_dict = defaultdict(dict)
        #
        self.thread = args.thread
        
        

    def split_fq(self,fq1,fq2,fq_outdir,sample,adapter_dict):

        out_notinfq1 = xopen(f"{fq_outdir}/{sample}_not_R1.fastq",'w')
        out_notinfq2 = xopen(f"{fq_outdir}/{sample}_not_R2.fastq",'w')

        out_fq1_dict = {adapter:xopen(f"{fq_outdir}/{sample}_{adapter}_R1.fastq",'w') for adapter in adapter_dict}
        out_fq2_dict = {adapter:xopen(f"{fq_outdir}/{sample}_{adapter}_R2.fastq",'w') for adapter in adapter_dict}
        
        with pysam.FastxFile(fq1) as fh1, pysam.FastxFile(fq2) as fh2:
            for record1,record2 in zip(fh1,fh2):
                header1, seq1, qual1 = record1.name, record1.sequence, record1.quality
                header2, seq2, qual2 = record2.name, record2.sequence, record2.quality
                
                adapter_distance = {adapter:hm_dis(seq2[:9],adapter_dict[adapter])  for adapter in adapter_dict}
                match_adapter = [k for k,v in adapter_distance.items() if v<=1]

                if len(match_adapter) == 0:
                    out_notinfq1.write(f'@{header1}\n{seq1[9:]}\n+\n{qual1[9:]}\n')
                    out_notinfq2.write(f'@{header2}\n{seq2[9:]}\n+\n{qual2[9:]}\n')    
                else:
                    out_fq1_dict[match_adapter[0]].write(f'@{header1}\n{seq1[9:]}\n+\n{qual1[9:]}\n')
                    out_fq2_dict[match_adapter[0]].write(f'@{header2}\n{seq2[9:]}\n+\n{qual2[9:]}\n')
        out_notinfq1.close()
        out_notinfq2.close()
        for adapter in adapter_dict:
            out_fq1_dict[adapter].close()
            out_fq2_dict[adapter].close()


    def run(self,params):
        fq1,fq2,outdir_fq,sample_name,adapter_dict = params
        self.split_fq(fq1,fq2,outdir_fq,sample_name,adapter_dict)


    def main(self):
        """
        split adapter
        """
        args_dict = self.args_dict
        
        param_list = []
        for sample in args_dict:
            sample_name = args_dict[sample]['sample_name']
            fa = args_dict[sample]['fa']
            fq1 = args_dict[sample]['fq1']
            fq2 = args_dict[sample]['fq2']
            outdir_fq = args_dict[sample]['fq_outdir']

            Path(outdir_fq).mkdir(parents=True,exist_ok=True)
            
            ## add fastq 
            if len(fq1.split(",")) > 1:
                tmp_fastq = Path(f'{self.outdir}/tmp_fastq')
                tmp_fastq.mkdir(parents=True,exist_ok=True)
                fq1_list = fq1.split(",")
                fq2_list = fq2.split(",")
                soft = "zcat" if fq1_list[0][-2:] == 'gz' else 'cat'
                fastq1 = " ".join(fq1_list)
                fastq2 = " ".join(fq2_list)
                cmd_r1 = f"{soft} {fastq1} |gzip > {tmp_fastq}/{sample_name}_R1.fastq.gz"
                cmd_r2 = f"{soft} {fastq2} |gzip > {tmp_fastq}/{sample_name}_R2.fastq.gz"
                debug_subprocess_call(cmd_r1)
                debug_subprocess_call(cmd_r2)
                fq1 = f'{tmp_fastq}/{sample_name}_R1.fastq.gz'
                fq2 = f'{tmp_fastq}/{sample_name}_R2.fastq.gz'
        

            adapter_dict = split_adapter(fa)
            
            for adapter in adapter_dict:
                self.units_dict[f"{sample_name}_{adapter}"] = {"sample":f"{sample_name}_{adapter}","unit":1,"platform":'ILLUMINA','fq1':f"{outdir_fq}/{sample_name}_{adapter}_R1.fastq",'fq2':f"{outdir_fq}/{sample_name}_{adapter}_R2.fastq"}
            
            param_list.append((fq1,fq2,outdir_fq,sample_name,adapter_dict))
        
        if self.thread > len(param_list):
            num_pool = len(param_list)
        else:
            num_pool = self.thread

        time1 = time.time()
        with multiprocessing.Pool(num_pool) as p:
            r=list(tqdm(p.map(self.run,param_list),total=len(param_list),desc='Split fastq '))
        p.close()
        p.join()
        time2 = time.time()
        used = timedelta(seconds=time2 - time1)
        print(f'use time {used}')
            
        self.units_config = pd.DataFrame.from_dict(self.units_dict).T
        self.units_config.to_csv(f"{self.config_path}/units.tsv",sep="\t",index=None)


@add_log
def split(args):
    runner =  Split_fastq(args)
    runner.main()
    
    
def get_opts_split(parser, sub_program=True):
    parser.add_argument('--mapfile',help='The necessary mapfile file location for splitting the original sequence.',type=str)
    parser.add_argument('--thread',help='Number of threads.', default=4,type=int)
    
    if sub_program:
        parser = s_common(parser)
    return parser


if __name__ == '__main__':
    unittest.main()