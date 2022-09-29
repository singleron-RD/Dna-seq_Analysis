import pysam
import pandas as pd
from pathlib import Path
import subprocess

from dna_seq_analysis.tools.common import *



GT_LIST_NOT = [(0, 0), (None, None)]

THREAD = 6



def get_only_child_dir(path):
    children = [child for child in path.iterdir() if child.is_dir()]
    assert (
        len(children) == 1
    ), "Invalid VEP cache directory, only a single entry is allowed, make sure that cache was created with the snakemake VEP cache wrapper"
    return children[0]


def fun(record):
    for _,value in record.samples.items():
        return value["GT"]


def repalce_html(html_file):
    """
    repalce html logo
    """
    JSAPI = CU_PATH.parents[2]/"js/jsapi"
    logo_file = CU_PATH.parents[2]/"templates/html/logo.html"
    with open(JSAPI, 'r') as file:
        replace_content = file.read()
        replace_content = "  "+'<script type="text/javascript">'+" "+replace_content+'</script>'
    with open(logo_file, 'r') as file:
        logo_content = file.read()
        
    sample_name = html_file.rsplit("/",1)[1].split(".")[0]
    
    output_html = Path(html_file).parents[0]/f'{sample_name}_sgr.html'
    
    with open(html_file, 'r') as file:
        fcontent = file.readlines()
        with open(output_html, 'w') as fp:
            for line in fcontent:
                if line.find('http://www.google.com/jsapi') != -1:
                    line = replace_content
                if line.find('<a href="http://www.ensembl.org/">') != -1:
                    line = logo_content
                if line.find('<a href="http://www.ensembl.org/vep">') != -1:
                    line = ""
                fp.write(line)
                
    cmd_clean_html = (f"rm {html_file}")
    subprocess.check_call(cmd_clean_html,shell=True)


class Annotate():
    """
    
    """
    def __init__(self,outdir):
        #self.sample,self.unit = wildcards
        self.outdir = Path(outdir)

        # input 
        self.plugins = f"{str(resource_dir)}/vep/plugins"
        self.cache = f"{str(resource_dir)}/vep/cache"
        
        # output
        self.calls = f"{str(self.outdir)}/filtered/all.vcf.gz"
        self.vcf2tsv_file = f"{str(self.outdir)}/tables/calls.tsv.gz"
        self.stats = f"{str(self.outdir)}/stats/all.stats.html"
        self.calls_output = f"{str(self.outdir)}/annotated/all.vcf.gz"

        annotated_dir = self.outdir/"annotated"
        annotated_dir.mkdir(parents=True,exist_ok=True)
        stats_dir = self.outdir/"stats"
        stats_dir.mkdir(parents=True,exist_ok=True)
        tables_dir = self.outdir/"tables"
        tables_dir.mkdir(parents=True,exist_ok=True)

        
    @add_log
    def annotation(self):
        """
        annotation all of sample
        """
        plugins = config["params"]["vep"]["plugins"]
        extra = config["params"]["vep"]["extra"]

        fork = "--fork {}".format(THREAD) if THREAD > 1 else ""
        load_plugins = " ".join(map("--plugin {}".format,plugins))

        if self.calls.endswith(".vcf.gz"):
            fmt = "z"
        elif self.calls.endswith(".bcf"):
            fmt = "b"
        else:
            fmt = "v"

        if self.cache:
            entrypath = get_only_child_dir(get_only_child_dir(Path(self.cache)))
            species = entrypath.parent.name
            release, build = entrypath.name.split("_")
            cache = (f"--offline --cache --dir_cache {self.cache} --cache_version {release} --species {species} --assembly {build}")

        cmd = (
            f"bcftools view {self.calls} | "
            f"vep {extra} {fork} "
            "--format vcf "
            f"--vcf {cache} "
            f"--dir_plugins {self.plugins} {load_plugins} "
            f"--output_file STDOUT --stats_file {self.stats} | "
            f"bcftools view -O{fmt} > {self.calls_output}"
        )
        debug_subprocess_call(cmd)
        repalce_html(self.stats)

    @add_log
    def vcf2tsv(self):
        """
        creat tables/calls.tsv.gz
        """
        cmd = (f"bcftools view --apply-filters PASS --output-type u {self.calls} | "
            "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
            f"gzip > {self.vcf2tsv_file}")
        debug_subprocess_call(cmd)
        
    def run(self):
        self.annotation()
        self.vcf2tsv()
        

class Split_vcf():
    """
    Each sample is split from the whole and its own report file is generated.
    """
    def __init__(self,outdir):
        self.outdir = Path(outdir)

        # input 
        self.plugins = f"{str(resource_dir)}/vep/plugins"
        self.cache = f"{str(resource_dir)}/vep/cache"
        self.calls = f"{str(self.outdir)}/annotated/all.vcf.gz"

    @add_log
    def split(self,sample_name):
        outdir = f"{str(self.outdir)}/split/{sample_name}"
        Path(outdir).mkdir(parents=True,exist_ok=True)
    
        cmd_split = (f'bcftools view -s {sample_name} {self.calls} -Oz -o {outdir}/{sample_name}_all.vcf.gz')
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

        plugins = config["params"]["vep"]["plugins"]
        extra = config["params"]["vep"]["extra"]

        fork = "--fork {}".format(THREAD) if THREAD > 1 else ""
        load_plugins = " ".join(map("--plugin {}".format,plugins))

        if self.calls.endswith(".vcf.gz"):
            fmt = "z"
        elif self.calls.endswith(".bcf"):
            fmt = "b"
        else:
            fmt = "v"

        if self.cache:
            entrypath = get_only_child_dir(get_only_child_dir(Path(self.cache)))
            species = entrypath.parent.name
            release, build = entrypath.name.split("_")
            cache = (f"--offline --cache --dir_cache {self.cache} --cache_version {release} --species {species} --assembly {build}")


        cmd_html = (
                    f"bcftools view {outdir}/{sample_name}.vcf.gz | "
                    f"vep {extra} {fork} --format vcf --vcf {cache} --dir_plugins {plugins} {load_plugins} --output_file STDOUT --stats_file {outdir}/{sample_name}.html | "
                    f"bcftools view -O{fmt} > {outdir}/{sample_name}_annotated.vcf.gz;"
                )
        subprocess.check_call(cmd_html,shell=True)
        
        html_file = f'{outdir}/{sample_name}.html'
        repalce_html(html_file)

        cmd_line_vcf2tsv = (
                    f"bcftools view --apply-filter PASS --output-type u {outdir}/{sample_name}_annotated.vcf.gz | "
                    f"rbt vcf-to-txt -g --fmt DP AD --info ANN | "
                    f"gzip > {outdir}/{sample_name}_calls.tsv.gz"
                    )
        subprocess.check_call(cmd_line_vcf2tsv,shell=True)


    def run_split_vcf(self):
        """
        split {sample}vcf.gz from all.vcf.gz 
        """
        calls_file = f"{str(self.outdir)}/tables/calls.tsv.gz"
        calls_df = pd.read_table(calls_file, header=[0, 1])
        samples = [name for name in calls_df.columns.levels[0] if name != "VARIANT"]
        for sample in samples:
            self.split(sample)

    @add_log
    def run(self):
        self.run_split_vcf()


def main():
    outdir = config['outdir']
    run_annotate = Annotate(outdir)
    run_annotate.run()
    run_split =  Split_vcf(outdir)
    run_split.run()


if __name__ == "__main__":
    main()