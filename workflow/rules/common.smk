import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
import pysam
from xopen import xopen
import re
from collections import defaultdict,OrderedDict
from tenacity import retry_if_result,retry
import pathlib
import subprocess



min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


##### Split_fastq ######
def split_adapter(fa):
    """
    creat adapter_dict
    {"adapter1":xxxxxxx}
    """
    adapter_dict ={}
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


def check_return_info(return_info):
    """
    check fun return
    """
    if return_info > 0:
        return True
    else:
        return False
    

class Split_fastq():
    """
    Split the original sequence and generate a Config file
    """
    def __init__(self,mapfile):
        self.mapfile = mapfile
        self.args_dict = parse_mapfile(self.mapfile)
    

    def run(self):
        """
        split adapter
        """
        args_dict = self.args_dict
        units_dict = defaultdict(dict)
        sample_dict = defaultdict(dict)
        for sample in args_dict:
            global I
            I = 0
            sample_name = args_dict[sample]['sample_name']
            fa = args_dict[sample]['fa']
            fq1 = args_dict[sample]['fq1']
            fq2 = args_dict[sample]['fq2']
            fq_outdir = args_dict[sample]['fq_outdir']
            config_outdir = "./config"
            mod = args_dict[sample]['mod']
        
            pathlib.Path(fq_outdir).mkdir(parents=True,exist_ok=True)
        
            adapter_dict = split_adapter(fa)
            adapter_dict_sorted = OrderedDict(sorted(adapter_dict.items(),key=lambda l: int(re.findall('\d+', l[0])[0])))
            
            index_len = len(adapter_dict_sorted)

            for adapter in adapter_dict:
                sample_dict[f"{sample_name}_{adapter}"] = {"sample":f"{sample_name}_{adapter}"}
                units_dict[f"{sample_name}_{adapter}"] = {"sample":f"{sample_name}_{adapter}","unit":1,"platform":'ILLUMINA','fq1':f"{fq_outdir}/{sample_name}_{adapter}_R1.fastq",'fq2':f"{fq_outdir}/{sample_name}_{adapter}_R2.fastq"}
              
            self.adapter_list(sample_name,fq1,fq2,fq_outdir,mod,adapter_dict_sorted,index_len)

        sample_config = pd.DataFrame.from_dict(sample_dict).T
        sample_config.to_csv(f"{config_outdir}/samples.tsv",sep="\t",index=None)
        units_config = pd.DataFrame.from_dict(units_dict).T
        units_config.to_csv(f"{config_outdir}/units.tsv",sep="\t",index=None)
            

    @retry(retry=retry_if_result(check_return_info))
    def adapter_list(self,sample_name,fq1,fq2,fq_outdir,mod,adapter_dict_sorted,index_len):
        global I
        if mod == 'fq1':
            fq=fq1
        if mod == 'fq2':
            fq=fq2

        frist_adapter = list(adapter_dict_sorted.items())[0][0]
        name_list = xopen(f"{fq_outdir}/{sample_name}_{frist_adapter}.lst",'w')

        fh =  pysam.FastxFile(fq)
        for record in fh:
            header, seq, = record.name, record.sequence
            if hm_dis(seq[:9],adapter_dict_sorted[frist_adapter]) > 1:
                continue
            name_list.write(f'{header}\n')
        name_list.close()
        fh.close()

        out_fq_R1 = f'{fq_outdir}/{sample_name}_{frist_adapter}_R1.fastq'
        out_fq_R2 = f'{fq_outdir}/{sample_name}_{frist_adapter}_R2.fastq'

        cmd_line = (
                    f'seqtk subseq {fq1} {fq_outdir}/{sample_name}_{frist_adapter}.lst > {out_fq_R1};'
                    f'seqtk subseq {fq2} {fq_outdir}/{sample_name}_{frist_adapter}.lst > {out_fq_R2}'
                    )
       
        subprocess.check_call(cmd_line,shell=True)
        print(f'the {I+1} finished')
        
        if I < index_len:
            adapter_dict_sorted.popitem(last=False)
            I += 1
            i = index_len - I 
            return i

###### Config file and sample sheets #####
configfile: "config/config.yaml"

if config["split_fq"]:
    mapfile = config['mapfile']
    runner =  Split_fastq(mapfile)
    runner.run()


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Helper functions #####

# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


def get_mapfile():
    """
    Get split fastq mapfile
    """
    if config["split_fq"]:
        # case 2: remove duplicates
        mapfile = "config/mapfile"
    return mapfile


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}-{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    )


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}
