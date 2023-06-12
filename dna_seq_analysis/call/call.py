from itertools import product
from pathlib import Path
import multiprocessing
from tqdm import tqdm
import unittest
import pandas as pd

from dna_seq_analysis.tools.common import *



class Calling():
    """
    calling vcf
    """
    def __init__(self,wildcards,outdir,chro,resource_dir,units,args):
        self.sample,_ = wildcards
        self.outdir = Path(outdir)
        self.chro = chro 
        self.resource_dir = resource_dir
        self.units = units
        self.args = args
        
        # input
        self.ref = f"{resource_dir}/genome.fasta"
        self.known = f"{resource_dir}/variation.noiupac.vcf.gz"
        self.known_idx = f"{resource_dir}/variation.noiupac.vcf.gz.tbi"

        # output
        called_dir = self.outdir/"05.called"
        called_dir.mkdir(parents=True,exist_ok=True)

        # log
        self.log_dir = self.outdir/"logs/called"
        self.log_dir.mkdir(parents=True,exist_ok=True)
        

    def call_variants(self):
        """
        """
        bams_list = get_sample_bams(self.units,self.outdir,self.sample)
        bams = ' '.join(list(map("-I {}".format, bams_list)))

        if self.known:
            known = "--dbsnp " + self.known
        
        if self.args.intervals:
            df = pd.read_csv(self.args.intervals,sep="\t",header=None)
            df.iloc[:,0] = df.iloc[:,0].astype('str')
            df[df.iloc[:,0] == str(self.chro)].to_csv(f"{self.outdir}/05.called/{self.chro}.regions.bed",sep="\t",header=None,index=None)

        extra = get_call_variants_params(self.chro,self.outdir,self.args.intervals,self.args.interval_padding)
        _log = self.log_dir/f"HaplotypeCaller-{self.sample}-{self.chro}.log"
        cmd = (
                f"gatk --java-options '-Xmx10g' HaplotypeCaller {extra} "
                f"-R {self.ref} {bams} "
                f"-ERC GVCF "
                f"-O {self.outdir}/05.called/{self.sample}.{self.chro}.g.vcf.gz {known} > {_log}"
                )
        if  not Path(f"{self.outdir}/05.called/{self.sample}.{self.chro}.g.vcf.gz").exists():
            try:
                debug_subprocess_call(cmd)
            except subprocess.CalledProcessError as e:
                print("Command failed with return code:", e.returncode)
                print(cmd)
                


    @add_log
    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.call_variants()


def run_call(params):
    """
    """
    wildcards,outdir,chro,resource_dir,units,args = params
    app_calling = Calling(wildcards,outdir,chro,resource_dir,units,args)
    app_calling.main()



class Combine():
    def __init__(self,outdir,chro,resource_dir,units):
        self.outdir = Path(outdir)
        self.chro = chro
        self.resource_dir = resource_dir
        self.units = units
        
        # input
        self.ref = f"{resource_dir}/genome.fasta"

        # output
        genotyped_dir = self.outdir/"06.genotyped"
        genotyped_dir.mkdir(parents=True,exist_ok=True)

        # log
        self.log_dir = self.outdir/"logs/called"


    def combine_calls(self):
        """
        """
        samples_set = set(self.units['sample'])
        _log = self.log_dir/f"CombineGVCFs-{self.chro}.log"
        gvcfs = [f"{str(self.outdir)}/05.called/{sample}.{self.chro}.g.vcf.gz" for sample in samples_set]
        gvcfs = ' '.join(list(map("-V {}".format, gvcfs)))
        cmd = (
            "gatk --java-options '-Xmx10g' CombineGVCFs "
            f"{gvcfs} "
            f"-R {self.ref} "
            f"-O {str(self.outdir)}/05.called/all.{self.chro}.g.vcf.gz > {_log}"
        )
        if  not Path(f"{str(self.outdir)}/05.called/all.{self.chro}.g.vcf.gz").exists():
            debug_subprocess_call(cmd)


    
    def genotype_variants(self):
        # Allow for either an input gvcf or GenomicsDB
        _log = self.log_dir/f"GenotypeGVCFs-{self.chro}.log"
        cmd=(
            f"gatk --java-options '-Xmx20g' GenotypeGVCFs  "
            f"-V {str(self.outdir)}/05.called/all.{self.chro}.g.vcf.gz "
            f"-R {self.ref} "
            f"-O {str(self.outdir)}/06.genotyped/all.{self.chro}.vcf.gz > {_log}"
        )
        if  not Path(f"{self.outdir}/06.genotyped/all.{self.chro}.vcf.gz").exists():
            debug_subprocess_call(cmd)


    @add_log
    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.combine_calls()
        self.genotype_variants()


def run_comine(params):
    """
    """
    outdir,chro,resource_dir,units = params
    app_combine = Combine(outdir,chro,resource_dir,units)
    app_combine.main()


@add_log
def merge_variants(outdir,chr_list):
    """
    """
    ## merge_variants
    outdir = Path(outdir)
    log_dir = outdir/"logs/called"
    _log = log_dir/f"MergeVcfs.log"
    inputs = " ".join(f"INPUT={str(outdir)}/06.genotyped/all.{chro}.vcf.gz" for chro in chr_list)
    cmd = (
        "picard MergeVcfs -Xmx4g "
        f"{inputs} "
        f"OUTPUT={str(outdir)}/06.genotyped/all.vcf.gz > {_log}"
    )
    debug_subprocess_call(cmd)


@add_log
def call(args):
    config_path = args.config_path
    threads = args.thread
    config = parse_config(config_path)
    outdir = config['outdir']
    resource_dir = config['genomedir']
    
    called_dir = Path(outdir)/"05.called"
    called_dir.mkdir(parents=True,exist_ok=True)
    bed =  Path(outdir)/"05.called/regions.bed"
    # bioawk -c fastx '{ print $name ":1-" length($seq) }
    cmd_bed = (
                "bioawk -c fastx '{ print $name }' "
                f"{resource_dir}/genome.fasta > {bed};"
                f"sed -i '/^[KGJ]/d' {bed}"
            )
    debug_subprocess_call(cmd_bed)
    
    units_file = f'{config_path}/units.tsv'
    units = get_units(units_file)
    
    # bed file
    chr_list = []
    with open(bed,'r') as fh:
        for chro in fh.readlines():
            chr_list.append(chro.strip())
    samples_set = set(units['sample'])
    
    # call
    call_param_list = []
    for wildcards in units.index:
        for chro in chr_list: 
            call_param_list.append((wildcards,outdir,chro,resource_dir,units,args))
    with multiprocessing.Pool(threads) as p:
        list(tqdm(p.map(run_call,call_param_list),total=len(call_param_list),desc='Calling '))
    p.close()
    p.join()

    # combine
    combine_param_list = []
    for chro in chr_list: 
        combine_param_list.append((outdir,chro,resource_dir,units))
    with multiprocessing.Pool(threads) as p:
        list(tqdm(p.map(run_comine,combine_param_list),total=len(combine_param_list),desc='Calling '))
    p.close()
    p.join()
    
    # merge all varians
    merge_variants(outdir,chr_list)

    # clean
    # factor Argument list too long
    for rm_file in [f"{str(outdir)}/05.called/{sample}.{chro}.g.vcf.gz" for sample,chro in product(samples_set,chr_list)]:
        filepath = Path(rm_file)
        filepath.unlink()
    for rm_file in [f"{str(outdir)}/06.genotyped/all.{chro}.vcf.gz" for chro in chr_list]:
        filepath = Path(rm_file)
        filepath.unlink()



def get_opts_call(parser, sub_program=True):
    parser.add_argument('--thread',help='Number of threads.', default=4,type=int)
    parser.add_argument('--intervals',help='One or more genomic intervals over which to operate.This argument may be specified 0 or more times.')
    parser.add_argument('--interval_padding',help='Amount of padding (in bp) to add to each interval you are including.',default=0,type=int)
    if sub_program:
        parser = s_common(parser)
    return parser


if __name__ == '__main__':
    unittest.main()