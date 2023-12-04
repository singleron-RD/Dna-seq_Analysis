import pysam
import pandas as pd
import numpy as np
from pathlib import Path
import os
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib
import unittest
from collections import Counter
import gseapy
from gseapy import barplot, dotplot
from tqdm import tqdm
import sys

from dna_seq_analysis.tools.common import *

matplotlib.use('Agg')

GT_LIST_NOT = [(0, 0), (None, None)]



def get_only_child_dir(path):
    children = [child for child in path.iterdir() if child.is_dir()]
    assert (
        len(children) == 1
    ), "Invalid VEP cache directory, only a single entry is allowed, make sure that cache was created with the ref module."
    return children[0]


def fun(record):
    for _,value in record.samples.items():
        return value["GT"]


def flatten(lst):
    for item in lst:
        if isinstance(item, list):
            yield from flatten(item)
        else:
            yield item


def repalce_html(html_file):
    """
    repalce html logo
    """
    JSAPI = CU_PATH.parents[1]/"templates/js/jsapi"
    logo_file = CU_PATH.parents[1]/"templates/html/logo.html"
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
    debug_subprocess_call(cmd_clean_html)


class Annotate():
    """
    
    """
    def __init__(self,outdir,resource_dir,args):
        #self.sample,self.unit = wildcards
        self.outdir = Path(outdir)
        self.args = args
        self.threads = args.annotation_thread
        # input 
        self.plugins = f"{resource_dir}/vep/plugins"
        self.cache = f"{resource_dir}/vep/cache"
        self.genome = f"{resource_dir}/genome.fasta"
        
        # output
        self.calls = f"{str(self.outdir)}/07.filtered/all.vcf.gz"
        self.vcf2tsv_file = f"{str(self.outdir)}/10.tables/calls.tsv.gz"
        self.stats = f"{str(self.outdir)}/09.stats/all.stats.html"
        self.calls_output = f"{str(self.outdir)}/08.annotated/all.vcf.gz"

        annotated_dir = self.outdir/"08.annotated"
        annotated_dir.mkdir(parents=True,exist_ok=True)
        stats_dir = self.outdir/"09.stats"
        stats_dir.mkdir(parents=True,exist_ok=True)
        tables_dir = self.outdir/"10.tables"
        tables_dir.mkdir(parents=True,exist_ok=True)

    def annotation(self):
        """
        annotation all of sample
        """
        plugins = self.args.vep_plugins_param
        extra = self.args.vep_param if self.args.vep_param else ''

        fork = "--fork {}".format(self.threads) if self.threads > 1 else ""
        load_plugins = f"--plugin {plugins},{self.plugins}/{plugins}_scores.txt"
        

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

        # use vep108
        cmd = (
            f"bcftools view {self.calls} | "
            f"vep {extra} {fork} "       
            "--format vcf "
            "--vcf --fields 'Allele,Consequence,SYMBOL,Gene,Feature,Feature_type,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons' "
            f"{cache} "
            f"--dir_plugins {self.plugins} {load_plugins} "
            "--force_overwrite "
            f"--hgvs --fasta {self.genome} "
            f"--output_file STDOUT --stats_file {self.stats} | "
            f"bcftools view -O{fmt} > {self.calls_output}"
        )
        debug_subprocess_call(cmd)
        repalce_html(self.stats)

    def vcf2tsv(self):
        """
        creat 10.tables/calls.tsv.gz
        """
        cmd = (f"bcftools view --apply-filters PASS --output-type u {self.calls} | "
            "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
            f"gzip > {self.vcf2tsv_file}")
        debug_subprocess_call(cmd)
    
    
    @add_log
    def run_all(self):
        self.annotation()
        self.vcf2tsv()
        

class Split_vcf():
    """
    Each sample is split from the whole and its own report file is generated.
    """
    def __init__(self,outdir,resource_dir,args):
        self.outdir = Path(outdir)
        self.args = args
        self.threads = args.annotation_thread
        self.species = args.species
        # input 
        self.plugins = f"{resource_dir}/vep/plugins"
        self.cache = f"{resource_dir}/vep/cache"
        self.calls = f"{str(self.outdir)}/08.annotated/all.vcf.gz"
        self.genome = f"{resource_dir}/genome.fasta"

    def split(self,sample_name):
        outdir = f"{str(self.outdir)}/11.split/{sample_name}"
        Path(outdir).mkdir(parents=True,exist_ok=True)
    
        cmd_split = (f'bcftools view -s {sample_name} {self.calls} -Oz -o {outdir}/{sample_name}_all.vcf.gz')
        debug_subprocess_call(cmd_split)

        with pysam.VariantFile(f"{outdir}/{sample_name}_all.vcf.gz") as vcf_in:
            with pysam.VariantFile(f"{outdir}/{sample_name}.vcf.gz",'w',header = vcf_in.header) as vcf_out:
                for record in vcf_in:
                    new_record = record.copy()
                    gt = fun(record)
                    if gt not in GT_LIST_NOT:
                        vcf_out.write(new_record)
        
        cmd_clean_vcf = (f'rm {outdir}/{sample_name}_all.vcf.gz')
        debug_subprocess_call(cmd_clean_vcf)

        plugins = self.args.vep_plugins_param
        extra = self.args.vep_param if self.args.vep_param else ''

        fork = "--fork {}".format(self.threads) if self.threads > 1 else ""
        load_plugins = f"--plugin {plugins},{self.plugins}/{plugins}_scores.txt"

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
                    f"vep {extra} {fork} "       
                    "--format vcf "
                    "--vcf --fields 'Allele,Consequence,SYMBOL,Gene,Feature,Feature_type,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons' "
                    f"{cache} "
                    f"--dir_plugins {self.plugins} {load_plugins} "
                    "--force_overwrite "
                    f"--hgvs --fasta {self.genome} "
                    f"--output_file STDOUT --stats_file {outdir}/{sample_name}.html | "
                    f"bcftools view -O{fmt} > {outdir}/{sample_name}_annotated.vcf.gz;"
                )
        debug_subprocess_call(cmd_html)
        
        html_file = f'{outdir}/{sample_name}.html'
        repalce_html(html_file)

        cmd_line_vcf2tsv = (
                    f"bcftools view --apply-filter PASS --output-type u {outdir}/{sample_name}_annotated.vcf.gz | "
                    f"rbt vcf-to-txt -g --fmt DP AD --info ANN | "
                    f"gzip > {outdir}/{sample_name}_calls.tsv.gz"
                    )
        debug_subprocess_call(cmd_line_vcf2tsv)
        
        # vcf2maf
        cmd_gzip = (f'gzip -d {outdir}/{sample_name}_annotated.vcf.gz')
        debug_subprocess_call(cmd_gzip)
        ncbi_build = self.args.build 
        cmd_vcf2maf = (f'vcf2maf.pl --input-vcf {outdir}/{sample_name}_annotated.vcf '
                       f'--output-maf {outdir}/{sample_name}_vep.maf '
                       f'--inhibit-vep --ref-fasta {self.genome} --ncbi-build {ncbi_build} '
                       f'--tumor-id {sample_name}'
                       )
        debug_subprocess_call(cmd_vcf2maf)
        
        # summary
        csq_type_all = []
        with pysam.VariantFile(f'{outdir}/{sample_name}_annotated.vcf','r') as vf:
            for record in vf:
                info = record.info
                csq = info['CSQ']
                csq_list = [x for x in csq]
                singel_all_type = [x.split("&") for x in [mix_x.split("|")[1] for mix_x in csq_list]]
                flat_singel_all_type = list((flatten(singel_all_type)))
                csq_type_all.extend(flat_singel_all_type)
        csq_dict = Counter(csq_type_all)
        df_csq = pd.DataFrame.from_dict(csq_dict,orient='index')
        df_csq.columns = ['Consequence type']
        df_csq.to_csv(f"{outdir}/{sample_name}_varitation_type.txt",sep="\t")
        
        
    def run_split_vcf(self):
        """
        split {sample}vcf.gz from all.vcf.gz 
        """
        calls_file = f"{str(self.outdir)}/10.tables/calls.tsv.gz"
        calls_df = pd.read_table(calls_file, header=[0, 1])
        samples = [name for name in calls_df.columns.levels[0] if name != "VARIANT"]
        
        p_list = []
        for sample in samples:
            p_list.append((sample))
        with multiprocessing.Pool(len(p_list)) as p:
            list(tqdm(p.imap(self.run_split,p_list),total=len(p_list),unit_scale = True,ncols = 70,file = sys.stdout,desc='Split vcf of each sample '))
        p.close()
        p.join()


    @add_log
    def run_split(self,param):
        sample = param
        self.split(sample)


def merge_maf(outdir):
    for sample in os.listdir(outdir):
        break
    cmd_extract_head = (f"grep 'Hugo_Symbol' {outdir}/{sample}/{sample}_vep.maf > {outdir}/tmp_head")
    cmd_extract_body = (f"cat {outdir}/*/*maf | grep -v '^#'| grep -v '^Hugo_Symbol' > {outdir}/tmp_all")
    cmd_merge = (f'cat {outdir}/tmp_head {outdir}/tmp_all > {outdir}/../08.annotated/vep_merge.maf')
    debug_subprocess_call(cmd_extract_head)
    debug_subprocess_call(cmd_extract_body)
    debug_subprocess_call(cmd_merge)
    Path(f'{outdir}/tmp_head').unlink()
    Path(f'{outdir}/tmp_all').unlink()
    

def maftools(outdir,species):
    cmd = (
        f'Rscript {str(CU_PATH.parents[1])}/tools/maftools.R '
        f'--outdir {outdir} '
        f'--species  {species} '
        f'--maf_file {outdir}/08.annotated/vep_merge.maf '
    )
    try:
        debug_subprocess_call(cmd)
    except subprocess.CalledProcessError as e:
        print("Command failed with return code:", e.returncode)


def go_kegg(outdir,species):
    if species == 'homo_sapiens':
        organism = 'Human'
        kegg_sets=['KEGG_2019_Human','KEGG_2021_Human']
        go_sets = ['GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021'] 
        en_dict = {'GO':go_sets,'KEGG':kegg_sets}
    elif species == 'mus_musculus':
        organism = 'mouse'
        kegg_sets=['KEGG_2019_Mouse']
        go_sets = ['GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021']
        en_dict = {'GO':go_sets,'KEGG':kegg_sets}
    
    df = pd.read_csv(f"{outdir}/08.annotated/mut_gene.txt",sep=",",header=0)
    sample_list = list(set(df['Tumor_Sample_Barcode']))
    for sample in sample_list:
        df_sub = df.query("Tumor_Sample_Barcode == @sample")
        gene_list = df_sub['Hugo_Symbol'].to_list()
        for enrichr_sampel in en_dict:
            try:
                enr = gseapy.enrichr(gene_list=gene_list,
                        gene_sets=en_dict[enrichr_sampel],
                        organism=organism,
                        outdir=f'{outdir}/11.split/{sample}',
                        cutoff=0.5
                        )
                # categorical scatterplot
                ax = dotplot(enr.results,
                            column="Adjusted P-value",
                            x='Gene_set',
                            size=10,
                            top_term=5,
                            figsize=(3,5),
                            title = "GO",
                            xticklabels_rot=45, 
                            show_ring=True,
                            marker='o',
                            cutoff=0.05,
                            ofname=f'{outdir}/11.split/{sample}/{sample}_{enrichr_sampel}_enrichr.pdf'
                            )      
                color_list=['darkred', 'darkblue','darkOrange']
                ax = barplot(enr.results,
                        column="Adjusted P-value",
                        group='Gene_set',
                        size=10,
                        top_term=5,
                        figsize=(3,5),
                        color=color_list, 
                        cutoff=0.05,
                        ofname=f'{outdir}/11.split/{sample}/{sample}_{enrichr_sampel}_enrichr.pdf'
                        )
            except ValueError:
                print(f'Warning: No enrich terms when cutoff = 0.05. of {sample}')


def var_type(x):
    if x[0] > x[1]:
        return 'deletion'
    elif (x[0] == x[1])&(x[0]==1):
        return 'SNV'
    elif (x[0] == x[1])&(x[0]!=1):
        return 'other'
    elif x[0] < x[1]:
        return 'insertion'


def plot_snv(split_dir):
    snv_dict = {}
    for sample in os.listdir(split_dir):
        cmd = (f"zcat {split_dir}/{sample}/{sample}_calls.tsv.gz|tail -n +3|cut -f 3,4 > {split_dir}/{sample}/{sample}.txt")
        debug_subprocess_call(cmd)
        df_sample = pd.read_csv(f"{split_dir}/{sample}/{sample}.txt",sep="\t",header=None)
        df_sample['len_tuble'] = df_sample.apply(lambda x: (len(x[0]),len(x[1])),axis=1)
        df_sample['var_type'] = df_sample['len_tuble'].apply(var_type)
        i = df_sample.var_type.value_counts().get('SNV',0)
        Path(f"{split_dir}/{sample}/{sample}.txt").unlink()
        snv_dict.update({sample:i})
    df = (
            pd.DataFrame.from_dict(snv_dict,orient='index')
            .sort_index()
            .reset_index()
        )
    df.columns = ['sample','total snvs']
    y_max = df['total snvs'].max()
    fig,ax = plt.subplots(figsize=(16,14))
    ax.bar(df['sample'],df['total snvs'],color='#0165fc',zorder=10)
    ax.set_ylim((0,y_max//5*6))
    ax.set_yticks(np.arange(0,y_max//5*7,step=y_max//5),[str(x) for x in np.arange(0,y_max//5*7,step=y_max//5)])
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.grid(which='major',axis='y',zorder=0)
    ax.tick_params(left=False,bottom=False)
    labels = ax.get_xticklabels()+ax.get_yticklabels()
    [label.set_fontname('serif') for label in labels]
    fontdict_lable={'size':17,
            'color':'k',
            'family':'serif'}
    ax.set_ylabel('Total SNVs',fontdict=fontdict_lable)
    # add title
    titlefountdict = {'size':20,
                      'color':'k',
                      'family':'serif'
                    }
    ax.set_title('High SNV detection rate',titlefountdict,pad=20)
    ax.set_xticklabels(labels=df['sample'].values.tolist(),rotation=15)
    fig.savefig(f"{split_dir}/total_snv.pdf",dpi=1000,bbox_inches='tight')
    fig.savefig(f"{split_dir}/total_snv.png",dpi=1000,bbox_inches='tight')



def annotation(args):
    config_path = args.config_path
    config = parse_config(config_path)
    outdir = config['outdir']
    resource_dir = config['genomedir']
    
    run_annotate = Annotate(outdir,resource_dir,args)
    run_annotate.run_all()    
    run_split_annotate =  Split_vcf(outdir,resource_dir,args)
    run_split_annotate.run_split_vcf()
    merge_maf(f'{outdir}/11.split')
    plot_snv(f'{outdir}/11.split')
    if args.species in ['homo_sapiens','mus_musculus']:
        maftools(outdir,args.species)
        go_kegg(outdir,args.species)
    else:
        print('Maftools and GO step is not running')


def get_opts_annotation(parser, sub_program=True):
    parser.add_argument('--annotation_thread',help='Number of threads in the annotation step.',type=int)
    parser.add_argument('--species',help="Ensembl species name.")
    parser.add_argument('--build',help="Genome build.")
    parser.add_argument('--vep_param',help="Additional parameters for the called software. Need to be enclosed in quotation marks.\
For example, `--{software}_param '--param1 value1 --param2 value2'`,`--vep_param '--sift b'`.")
    parser.add_argument('--vep_plugins_param',help="Add any plugin from https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html.Use named plugin.\
Multiple plugins can be used by supplying the --vep_plugins_param flag multiple times. For example, `--vep_plugins_param 'LoFtool'`.",default='LoFtool')
    
    if sub_program:
        parser = s_common(parser)
    return parser

if __name__ == "__main__":
    unittest.main()