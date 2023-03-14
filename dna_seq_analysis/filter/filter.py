from pathlib import Path
import unittest
from dna_seq_analysis.tools.common import *


HARD_DICT = {
    'snvs':"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
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


    def select_calls(self):
        for vartype in ['snvs','indels']:
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
        for vartype in ['snvs','indels']:
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
    

    def merge_calls(self):
        """
        """
        _log = self.log_dir/f"MergeVcfs.log"
        filtertype = self.args.filtertype
        inputs = " ".join(f"INPUT={str(self.outdir)}/07.filtered/all.{vartype}.{filtertype}.vcf.gz" for vartype in ['snvs','indels'])
        cmd = (
            "picard MergeVcfs -Xmx8g "
            f"{inputs} "
            f"OUTPUT={str(self.outdir)}/07.filtered/all.vcf.gz > {_log}"
        )
        debug_subprocess_call(cmd)

    @add_log
    def main(self):
        self.select_calls()
        self.hard_filter_calls()
        self.merge_calls()



@add_log
def filter(args):
    config_path = args.config_path
    config = parse_config(config_path)
    outdir = config['outdir']
    resource_dir = config['genomedir']
    
    app = Filtering(outdir,resource_dir,args)
    app.main()


def get_opts_filter(parser, sub_program=True):
    parser.add_argument('--filtertype',help='Types of filter.',choices=["recalibrated","hardfiltered"],default='recalibrated')
    
    if sub_program:
        parser = s_common(parser)

    return parser


if __name__ == '__main__':
    unittest.main()