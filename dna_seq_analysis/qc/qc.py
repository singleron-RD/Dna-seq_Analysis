from pathlib import Path
import multiprocessing
from tqdm import tqdm

from dna_seq_analysis.tools.common import *


class Fastqc():
    def __init__(self,fastq,outdir,wildcards):
        self.outdir = Path(outdir)
        self.sample,self.unit = wildcards
        self.fastq = fastq
        self.threads = 4

        # log
        log_dir = self.outdir/"logs"
        log_dir.mkdir(parents=True,exist_ok=True)
        self.fastqc_log = log_dir/'fastqc.log'

    def fastqc(self):
        """
        """
        # output
        fastqc_dir = self.outdir/"qc/fastqc"
        fastqc_dir.mkdir(parents=True,exist_ok=True)

        r1,r2 =  self.fastq['r1'], self.fastq['r2']
        cmd = (
        f"fastqc -t {self.threads} "
        f"--outdir {fastqc_dir} {r1} {r2} > {self.fastqc_log}"
    )
        debug_subprocess_call(cmd)


    def samtools_stats(self):
        """
        """
        # input 
        bam = self.outdir/f"recal/{self.sample}-{self.unit}.bam"
        # output
        samtools_stats_dir = self.outdir/"qc/samtools-stats"
        samtools_stats_dir.mkdir(parents=True,exist_ok=True)
        bam_txt = f"{str(samtools_stats_dir)}/{self.sample}-{self.unit}.txt"

        cmd = (f"samtools stats {bam} > {bam_txt} ")
        debug_subprocess_call(cmd)


    @add_log
    def main(self):
        self.fastqc()
        self.samtools_stats()


def run(params):
    fastq,outdir,wildcards = params
    app = Fastqc(fastq=fastq,outdir=outdir,wildcards=wildcards)
    app.main()
    


def main():
    outdir = config['outdir']
    param_list = []
    for wildcards in units.index:
        fastq = get_fastq(units,wildcards)
        param_list.append((fastq,outdir,wildcards))

    with multiprocessing.Pool(len(param_list)) as p:
        r=list(tqdm(p.map(run,param_list),total=len(param_list),desc='Fastqc '))
    p.close()
    p.join()

    fastqc_dir = Path(outdir)/"qc/fastqc"
    samtools_stats_dir = Path(outdir)/"qc/samtools-stats"
    dedup_dir = Path(outdir)/"qc/dedup"
    output_dir = Path(outdir)/"qc"

    multiqc_cmd = (
            "multiqc "
            "--force "
            f"-o {output_dir} -n multiqc.html "
            f"{fastqc_dir} {samtools_stats_dir} {dedup_dir}"
        )
    debug_subprocess_call(multiqc_cmd)


if __name__ == '__main__':
    main()

