import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import unittest

from dna_seq_analysis.tools.common import *


def stat(args):
    config_path = args.config_path
    config = parse_config(config_path)
    outdir = config['outdir']
    
    calls_tsv = f"{outdir}/10.tables/calls.tsv.gz"

    calls = pd.read_table(calls_tsv, header=[0, 1])
    samples = [name for name in calls.columns.levels[0] if name != "VARIANT"]
    sample_info = calls.loc[:, samples].stack([0, 1]).unstack().reset_index(1, drop=False)
    sample_info = sample_info.rename({"level_1": "sample"}, axis=1)

    sample_info = sample_info[sample_info["DP"] > 0]
    sample_info["freq"] = sample_info["AD"] / sample_info["DP"]
    sample_info.index = np.arange(sample_info.shape[0])

    # output
    plots_dir = Path(outdir)/"12.plots"
    plots_dir.mkdir(parents=True,exist_ok=True)

    freqs = f"{str(plots_dir)}/allele-freqs.svg"
    depths = f"{str(plots_dir)}/depths.svg"

    plt.figure()
    sns.stripplot(x="sample", y="freq", data=sample_info, jitter=True)
    plt.ylabel("allele frequency")
    plt.xticks(rotation="vertical")
    plt.savefig(freqs)

    plt.figure()
    sns.stripplot(x="sample", y="DP", data=sample_info, jitter=True)
    plt.ylabel("read depth")
    plt.xticks(rotation="vertical")
    plt.savefig(depths)
    
    
def get_opts_stat(parser, sub_program=True):
    if sub_program:
        parser = s_common(parser)
    return parser
    
    

if __name__ == "__main__":
    unittest.main()