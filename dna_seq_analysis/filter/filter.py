from pathlib import Path
import unittest
import pandas as pd


from dna_seq_analysis.tools.common import *


HARD_DICT = {
    'snps':"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
    'indels':"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
}


class Filtering():
    """
    """
    def __init__(self,outdir,resource_dir,args):
        self.outdir = Path(outdir)
        self.args = args
       
        # input
        self.ref = f"{resource_dir}/genome.fasta"
        self.known = f"{resource_dir}/variation.noiupac.vcf.gz"
        self.known_idx = f"{resource_dir}/variation.noiupac.vcf.gz.tbi"

        # output
        filtered_dir = self.outdir/"07.filtered"
        filtered_dir.mkdir(parents=True,exist_ok=True)

        # log
        self.log_dir = self.outdir/"logs/filtered"
        self.log_dir.mkdir(parents=True,exist_ok=True)
        
        if self.args.filtertype == 'recalibrated':
            if not self.args.species == 'homo_sapiens':
                raise ValueError("Only humans can choose `--filtertype recalibrated`,please change `--filtertype` parameter to hardfiltered.")     
        

    def select_calls(self):
        for vartype in ['snps','indels']:
            _log = self.log_dir/f"SelectVariants-{vartype}.log"
            extra = get_vartype_arg(vartype)
            cmd = (
                f"gatk --java-options '-Xmx4g' SelectVariants "
                f"-R {self.ref} "
                f"-V {str(self.outdir)}/06.genotyped/all.vcf.gz "
                f"{extra} "
                f"-O {str(self.outdir)}/07.filtered/all.{vartype}.vcf.gz > {_log}"
            )
            debug_subprocess_call(cmd)

    
    def hard_filter_calls(self):
        
        filtertype = self.args.filtertype
        for vartype in ['snps','indels']:
            _log = self.log_dir/f"VariantFiltration-{vartype}.log"
            filters = ''.join([
                "--filter-name {} --filter-expression '{}'".format(name, expr.replace("'", "\\'"))
                for name, expr in {"snv-hard-filter":HARD_DICT[vartype]}.items()
            ])
            cmd = (
                "gatk --java-options '-Xmx4g' VariantFiltration " 
                f"-R {self.ref} "
                f"-V {str(self.outdir)}/07.filtered/all.{vartype}.vcf.gz "
                f"{filters} "
                f"-O {str(self.outdir)}/07.filtered/all.{vartype}.{filtertype}.vcf.gz > {_log}"
            )
            debug_subprocess_call(cmd)
    
    
    def recalibrate_calls(self):
        # rename chr
        df = pd.DataFrame()
        chrom_list = list(range(1,23))+['X','Y','MT']
        df['old'] = chrom_list
        df['new'] = [f'chr{str(chm)}' for chm in chrom_list]
        df.iloc[24,1] = 'chrM'
        df.to_csv(f"{str(self.outdir)}/06.genotyped/chr_name_change_1.txt",sep="\t",index=None,header=None)
        df = df[['new','old']]
        df.to_csv(f"{str(self.outdir)}/06.genotyped/chr_name_change_2.txt",sep="\t",index=None,header=None)
        
        cmd_change_chrom_1 = (
                            f"bcftools annotate --rename-chrs {str(self.outdir)}/06.genotyped/chr_name_change_1.txt {str(self.outdir)}/06.genotyped/all.vcf.gz"
                            f"|bgzip -c > {str(self.outdir)}/06.genotyped/all_chr.vcf.gz"
        )
        debug_subprocess_call(cmd_change_chrom_1)
        cmd_index = (f"gatk IndexFeatureFile -I {str(self.outdir)}/06.genotyped/all_chr.vcf.gz")
        debug_subprocess_call(cmd_index)
        
        filtertype = self.args.filtertype 
        GATK_bundle = self.args.gatk_bundle_dir if self.args.gatk_bundle_dir else "/SGRNJ06/randd/public/wgs_ref/release-108/GATK_bundle/"
        REF = f"{GATK_bundle}/Homo_sapiens_assembly38.fasta"
        # SNP VQSR
        vartype = 'snps'
        mode = 'SNP'        
        _log = self.log_dir/f"VariantFiltration-{mode}.log"
        
        resources = f'--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {GATK_bundle}/hapmap_3.3.hg38.vcf.gz \
--resource:omini,known=false,training=true,truth=false,prior=12.0 {GATK_bundle}/1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 {GATK_bundle}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {GATK_bundle}/dbsnp_146.hg38.vcf.gz'

        annotation = "-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
        tranches = "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 "
        tranches_file = f"{str(self.outdir)}/07.filtered/all.{vartype}.tranches"
        tranches = "--tranches-file " + tranches_file
        cmd_snp = (
            f"gatk --java-options '-Xmx10G' VariantRecalibrator {resources} "
            f"-R {REF} -V {str(self.outdir)}/06.genotyped/all_chr.vcf.gz "
            f"-mode {mode} "
            f"--rscript-file {str(self.outdir)}/07.filtered/output.{vartype}.plots.R "
            f"-O {str(self.outdir)}/07.filtered/all.{vartype}.{filtertype}.vcf.gz "
            f"{tranches} {annotation} > {_log}"
        )
        
        apply_snp_vqsr_cmd = (f"gatk --java-options '-Xmx10G' ApplyVQSR "
                          f"-R {REF} -V {str(self.outdir)}/06.genotyped/all_chr.vcf.gz "
                          f"{tranches} "
                          f"--recal-file {str(self.outdir)}/07.filtered/all.{vartype}.{filtertype}.vcf.gz "
                          f"--truth-sensitivity-filter-level 99.0 "
                          "--max-gaussians 4 "
                          f"-mode {mode} "
                          f"-O {str(self.outdir)}/07.filtered/all.{vartype}.VQSR.vcf.gz")
        
        # Indel VQSR after SNP
        vartype = 'indels'
        mode = 'INDEL'        
        _log = self.log_dir/f"VariantFiltration-{mode}.log"
        
        resources = f'--resource:mills,known=true,training=true,truth=true,prior=12.0 {GATK_bundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {GATK_bundle}/dbsnp_146.hg38.vcf.gz'

        annotation = "-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
        tranches = "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 "
        tranches_file = f"{str(self.outdir)}/07.filtered/all.{vartype}.tranches"
        tranches = "--tranches-file " + tranches_file
        
        cmd_indel = (   f"gatk --java-options '-Xmx10G' VariantRecalibrator {resources} "
                        f"-R {REF} "
                        f"-V {str(self.outdir)}/07.filtered/all.snps.VQSR.vcf.gz "
                        f"-mode {mode} "
                        "--max-gaussians 4 "
                        f"--rscript-file {str(self.outdir)}/07.filtered/output.{vartype}.plots.R "
                        f"-O {str(self.outdir)}/07.filtered/all.{vartype}.VQSR.vcf.gz "
                        f"{tranches} {annotation} > {_log}"
        )
        apply_indel_vqsr_cmd = (f"gatk --java-options '-Xmx10G' ApplyVQSR "
                                f"-R {REF} "
                                f"-V {str(self.outdir)}/07.filtered/all.snps.VQSR.vcf.gz "
                                f"--truth-sensitivity-filter-level 99.0 "
                                f"{tranches} "
                                f"--recal-file {str(self.outdir)}/07.filtered/all.{vartype}.VQSR.vcf.gz "
                                f"-mode {mode} "
                                f"-O {str(self.outdir)}/07.filtered/all_chr.vcf.gz ")
        
        # run
        debug_subprocess_call(cmd_snp)
        debug_subprocess_call(apply_snp_vqsr_cmd)
        debug_subprocess_call(cmd_indel)
        debug_subprocess_call(apply_indel_vqsr_cmd)

        # rename
        cmd_change_chrom_2 = (
                            f"bcftools annotate --rename-chrs {str(self.outdir)}/06.genotyped/chr_name_change_2.txt {str(self.outdir)}/07.filtered/all_chr.vcf.gz"
                            f"|bgzip -c > {str(self.outdir)}/07.filtered/all.vcf.gz"
        )
        debug_subprocess_call(cmd_change_chrom_2)
        cmd_index = (f"gatk IndexFeatureFile -I {str(self.outdir)}/07.filtered/all.vcf.gz")
        debug_subprocess_call(cmd_index)

        # remove
        Path(f"{str(self.outdir)}/06.genotyped/chr_name_change_1.txt").unlink()
        Path(f"{str(self.outdir)}/06.genotyped/chr_name_change_2.txt").unlink()
        Path(f"{str(self.outdir)}/06.genotyped/all_chr.vcf.gz").unlink()
        Path(f"{str(self.outdir)}/06.genotyped/all_chr.vcf.gz.tbi").unlink()



    def merge_calls(self):
        """
        """
        _log = self.log_dir/"MergeVcfs.log"
        filtertype = self.args.filtertype
        inputs = " ".join(f"INPUT={str(self.outdir)}/07.filtered/all.{vartype}.{filtertype}.vcf.gz" for vartype in ['snps','indels'])
      
        cmd = (
            "picard MergeVcfs -Xmx8g "
            f"{inputs} "
            f"OUTPUT={str(self.outdir)}/07.filtered/all.vcf.gz > {_log}"
        )
        debug_subprocess_call(cmd)

    @add_log
    def main(self):
        if self.args.filtertype == "hardfiltered":
            self.select_calls()
            self.hard_filter_calls()
            self.merge_calls()
        else:
            self.recalibrate_calls()      
            


def filter(args):
    config_path = args.config_path
    config = parse_config(config_path)
    outdir = config['outdir']
    resource_dir = config['genomedir']
    
    app = Filtering(outdir,resource_dir,args)
    app.main()


def get_opts_filter(parser, sub_program=True):
    parser.add_argument('--filtertype',help='Types of filter.',choices=["recalibrated","hardfiltered"],default='hardfiltered')
    parser.add_argument('--species',help="Ensembl species name.")
    parser.add_argument('--gatk_bundle_dir',help="GATK bundle dir path.")
    
    
    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == '__main__':
    unittest.main()