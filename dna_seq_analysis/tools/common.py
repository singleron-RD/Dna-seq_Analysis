import subprocess
import logging
import time
from datetime import timedelta
from functools import wraps
import yaml
import pandas as pd
from pathlib import Path
import argparse


###### date #####

CU_PATH = Path(__file__).absolute()

parser = argparse.ArgumentParser(description='dna-seq_analysis')
parser.add_argument('--config_path',help='The position where the config.yaml is located ',required=True)
args = parser.parse_args()
config_path = args.config_path



config_path = config_path
configfile = Path(config_path)/"config.yaml"
with open(configfile,"r") as test_file:
    config = yaml.load(test_file,Loader=yaml.FullLoader)

if Path(config["units"]).exists():
    units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)

resource_dir = Path(config['resource'])

###### function #####
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


def debug_subprocess_call(cmd):
    '''
    debug subprocess call
    '''
    if cmd.find('2>&1') == -1:
        cmd += ' 2>&1 '
    subprocess.check_call(cmd, shell=True)


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


def get_fastq(units,wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[wildcards, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default



def get_read_group(units,sample,unit):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=sample,
        platform=units.loc[(sample, unit), "platform"],
    )



def get_recal_input(sample,unit,bai=False):
    # case 1: no duplicate removal
    f = f"mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = f"dedup/{sample}-{unit}.bam"
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


def get_sample_bams(units,dirs,sample):
    """Get all aligned reads of given sample."""

    unit=units.loc[sample]['unit']
    bams = [f"{str(dirs)}/recal/{sample}-{num}.bam" for num in unit]
    return bams



def get_call_variants_params(contig):
    regions = f"results/called/{contig}.regions.bed" if config["processing"].get("restrict-regions") else ""
    return (
        get_regions_param(
            regions=regions, default="--intervals {}".format(contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_vartype_arg(vartype):
    return "--select-type-to-include {}".format("SNP" if vartype == "snvs" else "INDEL")


def get_filter(vartype):
    return {"snv-hard-filter": config["filtering"]["hard"][vartype]}