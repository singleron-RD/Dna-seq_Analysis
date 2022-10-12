from itertools import product
from pathlib import Path
import multiprocessing
from tqdm import tqdm

from dna_seq_analysis.tools.common import *



class Calling():
    """
    calling vcf
    """
    def __init__(self,wildcards,outdir,chro):
        self.sample,_ = wildcards
        self.outdir = Path(outdir)
        self.chro = chro 
        
        # input
        self.ref = str(resource_dir/"genome.fasta")
        self.known = str(resource_dir/"variation.noiupac.vcf.gz")
        self.known_idx = str(resource_dir/"variation.noiupac.vcf.gz.tbi")

        # output
        called_dir = self.outdir/"called"
        called_dir.mkdir(parents=True,exist_ok=True)

        genotyped_dir = self.outdir/"genotyped"
        genotyped_dir.mkdir(parents=True,exist_ok=True)

        # log
        self.log_dir = self.outdir/"logs/called"
        self.log_dir.mkdir(parents=True,exist_ok=True)
        

    def call_variants(self):
        """
        """
        bams_list = get_sample_bams(units,self.outdir,self.sample)
        bams = ' '.join(list(map("-I {}".format, bams_list)))

        if self.known:
            known = "--dbsnp " + self.known

        
        extra = get_call_variants_params(self.chro)
        _log = self.log_dir/f"HaplotypeCaller-{self.sample}-{self.chro}.log"
        cmd = (
                f"gatk --java-options '-Xmx10g' HaplotypeCaller {extra} "
                f"-R {self.ref} {bams} "
                f"-ERC GVCF "
                f"-O {self.outdir}/called/{self.sample}.{self.chro}.g.vcf.gz {known} > {_log}"
                )
        if  not Path(f"{self.outdir}/called/{self.sample}.{self.chro}.g.vcf.gz").exists():
            debug_subprocess_call(cmd)


    @add_log
    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.call_variants()


def run_call(params):
    """
    """
    wildcards,outdir,chro = params
    app_calling = Calling(wildcards,outdir,chro)
    app_calling.main()


class Combine():
    def __init__(self,outdir,chro):
        self.outdir = Path(outdir)
        self.chro = chro
        
        # input
        self.ref = str(resource_dir/"genome.fasta")

        # output
        genotyped_dir = self.outdir/"genotyped"
        genotyped_dir.mkdir(parents=True,exist_ok=True)

        # log
        self.log_dir = self.outdir/"logs/called"


    @add_log
    def combine_calls(self):
        """
        """
        samples_set = set(units['sample'])
        _log = self.log_dir/f"CombineGVCFs-{self.chro}.log"
        gvcfs = [f"{str(self.outdir)}/called/{sample}.{self.chro}.g.vcf.gz" for sample in samples_set]
        gvcfs = ' '.join(list(map("-V {}".format, gvcfs)))
        cmd = (
            "gatk --java-options '-Xmx10g' CombineGVCFs "
            f"{gvcfs} "
            f"-R {self.ref} "
            f"-O {str(self.outdir)}/called/all.{self.chro}.g.vcf.gz > {_log}"
        )
        if  not Path(f"{str(self.outdir)}/called/all.{self.chro}.g.vcf.gz").exists():
            debug_subprocess_call(cmd)


    @add_log
    def genotype_variants(self):
        """
        """
        extra = config["params"]["gatk"]["GenotypeGVCFs"]
        # Allow for either an input gvcf or GenomicsDB
        _log = self.log_dir/f"GenotypeGVCFs-{self.chro}.log"
        cmd=(
            f"gatk --java-options '-Xmx20g' GenotypeGVCFs {extra} "
            f"-V {str(self.outdir)}/called/all.{self.chro}.g.vcf.gz "
            f"-R {self.ref} "
            f"-O {str(self.outdir)}/genotyped/all.{self.chro}.vcf.gz > {_log}"
        )
        if  not Path(f"{self.outdir}/genotyped/all.{self.chro}.vcf.gz").exists():
            debug_subprocess_call(cmd)


    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.combine_calls()
        self.genotype_variants()


def run_comine(params):
    """
    """
    outdir,chro = params
    app_combine = Combine(outdir,chro)
    app_combine.main()


@add_log
def merge_variants(outdir,chr_list):
    """
    """
    ## merge_variants
    outdir = Path(outdir)
    log_dir = outdir/"logs/called"
    _log = log_dir/f"MergeVcfs.log"
    inputs = " ".join(f"INPUT={str(outdir)}/genotyped/all.{chro}.vcf.gz" for chro in chr_list)
    cmd = (
        "picard MergeVcfs -Xmx4g "
        f"{inputs} "
        f"OUTPUT={str(outdir)}/genotyped/all.vcf.gz > {_log}"
    )
    debug_subprocess_call(cmd)


def main():
    outdir = config['outdir']
    called_dir = Path(outdir)/"called"
    called_dir.mkdir(parents=True,exist_ok=True)
    bed =  Path(outdir)/"called/regions.bed"
    cmd_bed = (
                "bioawk -c fastx '{ print $name }' "
                f"{str(resource_dir)}/genome.fasta > {bed}"
            )
    debug_subprocess_call(cmd_bed)

    # thread
    thread = 10
    
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
            call_param_list.append((wildcards,outdir,chro))
    with multiprocessing.Pool(thread) as p:
        r=list(tqdm(p.map(run_call,call_param_list),total=len(call_param_list),desc='Calling '))
    p.close()
    p.join()

    # combine
    combine_param_list = []
    for chro in chr_list: 
        combine_param_list.append((outdir,chro))
    with multiprocessing.Pool(thread) as p:
        r=list(tqdm(p.map(run_comine,combine_param_list),total=len(combine_param_list),desc='Calling '))
    p.close()
    p.join()
    
    # merge all varians
    merge_variants(outdir,chr_list)

    # clean
    # factor Argument list too long
    for rm_file in [f"{str(outdir)}/called/{sample}.{chro}.g.vcf.gz" for sample,chro in product(samples_set,chr_list)]:
        filepath = Path(rm_file)
        filepath.unlink()
    for rm_file in [f"{str(outdir)}/genotyped/all.{chro}.vcf.gz" for chro in chr_list]:
        filepath = Path(rm_file)
        filepath.unlink()



if __name__ == '__main__':
    main()