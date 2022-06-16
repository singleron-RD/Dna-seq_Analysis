import pysam
import pandas as pd
import os
import argparse
from pathlib import Path
import subprocess
import logging
import time
from datetime import timedelta
from functools import wraps



GT_LIST_NOT = [(0, 0), (None, None)]
current_dir = os.path.dirname(os.path.abspath(__file__))
workdir = current_dir.rsplit("/",1)[0]
FLODER = Path(workdir)



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


def fun(record):
    for _,value in record.samples.items():
        return value["GT"]



class Split_vcf():
    """
    Each sample is split from the whole and its own report file is generated.
    """
    def __init__(self,qc):
        self.qc = qc
        

    @add_log
    def split(self,sample_name):
        all_vcf = FLODER/"results/filtered/all.vcf.gz"
        outdir = FLODER/f"split/{sample_name}"
        Path(outdir).mkdir(parents=True,exist_ok=True)
    
        cmd_split = (f'bcftools view -s {sample_name} {all_vcf} -Oz -o {outdir}/{sample_name}_all.vcf.gz')
        subprocess.check_call(cmd_split,shell=True)

        with pysam.VariantFile(f"{outdir}/{sample_name}_all.vcf.gz") as vcf_in:
            with pysam.VariantFile(f"{outdir}/{sample_name}.vcf.gz",'w',header = vcf_in.header) as vcf_out:
                for record in vcf_in:
                    new_record = record.copy()
                    gt = fun(record)
                    if gt not in GT_LIST_NOT:
                        vcf_out.write(new_record)
        
        cmd_clean_vcf = (f'rm {outdir}/{sample_name}_all.vcf.gz')
        subprocess.check_call(cmd_clean_vcf,shell=True)

        cmd_html = (
                    f"bcftools view {outdir}/{sample_name}.vcf.gz |vep --fork 4 --format vcf --vcf --offline --cache --dir_cache /SGRNJ06/randd/public/wgs_ref/homo_sapiens/vep/cache --cache_version 98 --species homo_sapiens --assembly GRCh38 --dir_plugins /SGRNJ06/randd/public/wgs_ref/homo_sapiens/vep/plugins --plugin L --plugin o --plugin F --plugin t --plugin o --plugin o --plugin l --output_file STDOUT --stats_file {outdir}/{sample_name}.html| bcftools view -Oz > {outdir}/{sample_name}_annotated.vcf.gz;"
                )
        subprocess.check_call(cmd_html,shell=True)

        JSAPI = FLODER/"js/jsapi"
        with open(JSAPI, 'r') as file:
            replace_content = file.read()
            replace_content = "  "+'<script type="text/javascript">'+" "+replace_content+'</script>'

        html_file = f'{outdir}/{sample_name}.html'
        with open(html_file, 'r') as file:
            fcontent = file.readlines()
            output_html = f"{outdir}/{sample_name}_read.html"
            with open(output_html, 'w') as fp:
                for line in fcontent:
                    if line.find('http://www.google.com/jsapi') != -1:
                        line = replace_content
                    fp.write(line)
                    
        cmd_clean_html = (f"rm {html_file}")
        subprocess.check_call(cmd_clean_html,shell=True)

        cmd_line_vcf2tsv = (
                    f"bcftools view --apply-filter PASS --output-type u {outdir}/{sample_name}_annotated.vcf.gz | "
                    f"rbt vcf-to-txt -g --fmt DP AD --info ANN | "
                    f"gzip > {outdir}/{sample_name}_calls.tsv.gz"
                    )
        subprocess.check_call(cmd_line_vcf2tsv,shell=True)


    @add_log
    def run_split_vcf(self):
        """
        split {sample}vcf.gz from all.vcf.gz 
        """
        calls_file = FLODER/"results/tables/calls.tsv.gz"
        calls_df = pd.read_table(calls_file, header=[0, 1])
        samples = [name for name in calls_df.columns.levels[0] if name != "VARIANT"]
        for sample in samples:
            self.split(sample)
        if self.qc == "yes":
            qc_sh = current_dir+"/collectwgs.sh"
            recal_dir = FLODER/"results/recal"
            cmd_cp = (f'cp {qc_sh} {recal_dir}')
            subprocess.check_call(cmd_cp,shell=True)
            os.chdir(FLODER/"results/recal")
            cmd_qc = (f'sh collectwgs.sh')
            subprocess.check_call(cmd_qc,shell=True)
        

    @add_log
    def run(self):
        self.run_split_vcf()
                 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--qc',help='Whether it is necessary to run quality control script to generate coverage and other indicators.',required=True,choices=['yes', 'no'])
    args = parser.parse_args()
    qc = args.qc
    runner =  Split_vcf(qc)
    runner.run()
    