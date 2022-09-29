from pathlib import Path
from dna_seq_analysis.tools.common import *



class Filtering():
    """
    """
    def __init__(self,outdir):
        self.outdir = Path(outdir)
        
        # input
        self.ref = str(resource_dir/"genome.fasta")
        self.known = str(resource_dir/"variation.noiupac.vcf.gz")
        self.known_idx = str(resource_dir/"variation.noiupac.vcf.gz.tbi")

        # output
        filtered_dir = self.outdir/"filtered"
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
                f"-V {str(self.outdir)}/genotyped/all.vcf.gz "
                f"{extra} "
                f"-O {str(self.outdir)}/filtered/all.{vartype}.vcf.gz > {_log}"
            )
            debug_subprocess_call(cmd)

    
    def hard_filter_calls(self):
        
        for vartype in ['snvs','indels']:
            _log = self.log_dir/f"VariantFiltration-{vartype}.log"
            filters = ''.join([
                "--filter-name {} --filter-expression '{}'".format(name, expr.replace("'", "\\'"))
                for name, expr in get_filter(vartype).items()
            ])
            cmd = (
                "gatk --java-options '-Xmx4g' VariantFiltration " 
                f"-R {self.ref} "
                f"-V {str(self.outdir)}/filtered/all.{vartype}.vcf.gz "
                f"{filters} "
                f"-O {str(self.outdir)}/filtered/all.{vartype}.hardfiltered.vcf.gz > {_log}"
            )
            debug_subprocess_call(cmd)
    

    def recalibrate_calls(self):
        """
        gatk VariantRecalibrator
        gatk ApplyVQSR


        extra = snakemake.params.get("extra", "")
        java_opts = get_java_opts(snakemake)

        def fmt_res(resname, resparams):
            fmt_bool = lambda b: str(b).lower()
            try:
                f = snakemake.input.get(resname)
            except KeyError:
                raise RuntimeError(
                    "There must be a named input file for every resource (missing: {})".format(
                        resname
                    )
                )
            return "{},known={},training={},truth={},prior={} {}".format(
                resname,
                fmt_bool(resparams["known"]),
                fmt_bool(resparams["training"]),
                fmt_bool(resparams["truth"]),
                resparams["prior"],
                f,
            )

        resources = [
            "--resource:{}".format(fmt_res(resname, resparams))
            for resname, resparams in snakemake.params["resources"].items()
        ]
        annotation = list(map("-an {}".format, snakemake.params.annotation))
        tranches = ""
        if snakemake.output.tranches:
            tranches = "--tranches-file " + snakemake.output.tranches

        log = snakemake.log_fmt_shell(stdout=True, stderr=True)
        shell(
            "gatk --java-options '{java_opts}' VariantRecalibrator {extra} {resources} "
            "-R {snakemake.input.ref} -V {snakemake.input.vcf} "
            "-mode {snakemake.params.mode} "
            "--output {snakemake.output.vcf} "
            "{tranches} {annotation} {log}"
        )
        """
        pass


    def merge_calls(self):
        """
        """
        _log = self.log_dir/f"MergeVcfs.log"
        filtertype = "recalibrated" if config["filtering"]["vqsr"] else "hardfiltered"
        inputs = " ".join(f"INPUT={str(self.outdir)}/filtered/all.{vartype}.{filtertype}.vcf.gz" for vartype in ['snvs','indels'])
        cmd = (
            "picard MergeVcfs -Xmx8g "
            f"{inputs} "
            f"OUTPUT={str(self.outdir)}/filtered/all.vcf.gz > {_log}"
        )
        debug_subprocess_call(cmd)

    @add_log
    def main(self):
        """
        """
        self.select_calls()
        self.hard_filter_calls()
        self.recalibrate_calls()
        self.merge_calls()



def main():
    outdir = config['outdir']
    app = Filtering(outdir)
    app.main()


if __name__ == '__main__':
    main()