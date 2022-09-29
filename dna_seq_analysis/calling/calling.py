from itertools import product
from pathlib import Path
import multiprocessing
from tqdm import tqdm

from dna_seq_analysis.tools.common import *



class Calling():
    """
    calling vcf
    """
    def __init__(self,wildcards,outdir,bed):
        self.sample,_ = wildcards
        self.outdir = Path(outdir)
        self.bed = bed
        
        # input
        self.ref = str(resource_dir/"genome.fasta")
        self.known = str(resource_dir/"variation.noiupac.vcf.gz")
        self.known_idx = str(resource_dir/"variation.noiupac.vcf.gz.tbi")

        # output
        called_dir = self.outdir/"called"
        called_dir.mkdir(parents=True,exist_ok=True)

        genotyped_dir = self.outdir/"genotyped"
        genotyped_dir.mkdir(parents=True,exist_ok=True)

        # bed file
        self.chr_list = []
        with open(self.bed,'r') as fh:
            for chro in fh.readlines():
                self.chr_list.append(chro.strip())

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

        for chro in self.chr_list:
            extra = get_call_variants_params(chro)
            _log = self.log_dir/f"HaplotypeCaller-{self.sample}-{chro}.log"
            cmd = (
                    f"gatk --java-options '-Xmx10g' HaplotypeCaller {extra} "
                    f"-R {self.ref} {bams} "
                    f"-ERC GVCF "
                    f"-O {self.outdir}/called/{self.sample}.{chro}.g.vcf.gz {known} > {_log}"
                    )
            if  not Path(f"{self.outdir}/called/{self.sample}.{chro}.g.vcf.gz").exists():
                debug_subprocess_call(cmd)


    @add_log
    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.call_variants()


class Combine():
    def __init__(self,outdir,bed):
        self.outdir = Path(outdir)
        self.bed = bed
        
        # input
        self.ref = str(resource_dir/"genome.fasta")


        # output
        genotyped_dir = self.outdir/"genotyped"
        genotyped_dir.mkdir(parents=True,exist_ok=True)

        # bed file
        self.chr_list = []
        with open(self.bed,'r') as fh:
            for chro in fh.readlines():
                self.chr_list.append(chro.strip())

        # log
        self.log_dir = self.outdir/"logs/called"


    @add_log
    def combine_calls(self):
        """
        """
        samples_set = set(units['sample'])
        for chro in self.chr_list:
            _log = self.log_dir/f"CombineGVCFs-{chro}.log"
            gvcfs = [f"{str(self.outdir)}/called/{sample}.{chro}.g.vcf.gz" for sample in samples_set]
            gvcfs = ' '.join(list(map("-V {}".format, gvcfs)))
            cmd = (
                "gatk --java-options '-Xmx10g' CombineGVCFs "
                f"{gvcfs} "
                f"-R {self.ref} "
                f"-O {str(self.outdir)}/called/all.{chro}.g.vcf.gz > {_log}"
            )
            debug_subprocess_call(cmd)


    @add_log
    def genotype_variants(self):
        """
        """
        extra = config["params"]["gatk"]["GenotypeGVCFs"]
        # Allow for either an input gvcf or GenomicsDB
        for chro in self.chr_list:
            _log = self.log_dir/f"GenotypeGVCFs-{chro}.log"
            cmd=(
                f"gatk --java-options '-Xmx20g' GenotypeGVCFs {extra} "
                f"-V {str(self.outdir)}/called/all.{chro}.g.vcf.gz "
                f"-R {self.ref} "
                f"-O {str(self.outdir)}/genotyped/all.{chro}.vcf.gz > {_log}"
            )
            if  not Path(f"{self.outdir}/genotyped/all.{chro}.vcf.gz").exists():
                debug_subprocess_call(cmd)


    @add_log
    def merge_variants(self):
        """
        """
        ## merge_variants
        _log = self.log_dir/f"MergeVcfs.log"
        inputs = " ".join(f"INPUT={str(self.outdir)}/genotyped/all.{chro}.vcf.gz" for chro in self.chr_list)
        cmd = (
            "picard MergeVcfs -Xmx4g "
            f"{inputs} "
            f"OUTPUT={str(self.outdir)}/genotyped/all.vcf.gz > {_log}"
        )
        debug_subprocess_call(cmd)


    def main(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.combine_calls()
        self.genotype_variants()
        self.merge_variants()



def run(params):
    """
    """
    wildcards,outdir,bed = params
    app_calling = Calling(wildcards,outdir,bed)
    app_calling.main()




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

    param_list = []
    for wildcards in units.index: 
        param_list.append((wildcards,outdir,bed))

    with multiprocessing.Pool(len(param_list)) as p:
        r=list(tqdm(p.map(run,param_list),total=len(param_list),desc='Calling '))
    p.close()
    p.join()

    app_combine = Combine(outdir,bed)
    app_combine.main()

    # clean
    chr_list = []
    with open(bed,'r') as fh:
        for chro in fh.readlines():
            chr_list.append(chro.strip())
    samples_set = set(units['sample'])
    called_gvcfs = " ".join([f"{str(outdir)}/called/{sample}.{chro}.g.vcf.gz" for sample,chro in product(samples_set,chr_list)])
    genotyped_gvcfs = " ".join(f"{str(outdir)}/genotyped/all.{chro}.vcf.gz" for chro in chr_list)
    cmd_clean = (f"rm {called_gvcfs} {genotyped_gvcfs}")
    debug_subprocess_call(cmd_clean)



if __name__ == '__main__':
    main()