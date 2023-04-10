import subprocess
import logging
import time
from datetime import timedelta
from functools import wraps
import yaml
import pandas as pd
from pathlib import Path
import pysam
import pathlib
import os
import importlib

from dna_seq_analysis.__init__ import ROOT_PATH


CU_PATH = Path(__file__).absolute()


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


def s_common(parser):
    """subparser common arguments
    """
    parser.add_argument('--config_path',help='The position where the config.yaml is located ',required=True)
    return parser



def parse_config(config_path):
    configfile = Path(config_path)/"config.yaml"
    with open(configfile,"r") as arg_fh:
        config = yaml.load(arg_fh,Loader=yaml.FullLoader)
    return config


def get_units(units_file):
    if Path(units_file).exists():
        units = pd.read_table(units_file, dtype=str).set_index(["sample", "unit"], drop=False)
    return units


def find_step_module(step):
    file_path_dict = {
        'step': f'{ROOT_PATH}/{step}/{step}.py',
    }
    if os.path.exists(file_path_dict['step']):
        step_module = importlib.import_module(f"dna_seq_analysis.{step}.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {step}.{step}")

    return step_module



class Step:
    """
    Step class
    """

    def __init__(self, args, display_title=None):
        '''
        display_title controls the section title in HTML report
        '''
        print(f'Args: {args}')
        self.args = args
        self.display_title = display_title
        


def debug_subprocess_call(cmd):
    '''
    debug subprocess call
    '''
    if cmd.find('2>&1') == -1:
        cmd += ' 2>&1 '
    subprocess.check_call(cmd, shell=True)


def get_fastq(units,wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[wildcards, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_read_group(units,sample,unit):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=sample,
        platform=units.loc[(sample, unit), "platform"],
    )


def get_recal_input(sample,unit,remove_duplicates,bai=False,restrict_regions=None):
    # case 1: no duplicate removal
    f = f"02.mapped/{sample}-{unit}.sorted.bam"
    if remove_duplicates:
        # case 2: remove duplicates
        f = f"03.dedup/{sample}-{unit}.bam"
    if bai:
        if restrict_regions:
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
    bams = [f"{str(dirs)}/04.recal/{sample}-{num}.bam" for num in unit]
    return bams


def get_bam_file(file_path, pattern="*.bam"):
    all_file = []
    files = pathlib.Path(file_path).rglob(pattern)
    for file in files:
        if pathlib.Path.is_file(file):
            all_file.append(file)
    return all_file


def get_regions_param(regions,padding,default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(contig,dirs,restrict_regions,restrict_padding):
    regions = f"{dirs}/05.called/{contig}.regions.bed" if restrict_regions else ""
    return (
        get_regions_param(
            regions=regions,padding=restrict_padding, default="--intervals {}".format(contig)
        )
    )


def get_vartype_arg(vartype):
    return "--select-type-to-include {}".format("SNP" if vartype == "snvs" else "INDEL")