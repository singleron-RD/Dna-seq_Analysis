import os
from pathlib import Path
import multiprocessing
from tqdm import tqdm
from collections import defaultdict
import pysam

from dna_seq_analysis.tools.common import *



class Map_reads():
    """
    bwa mapping
    """
    def __init__(self,wildcards,threads,outdir):
        self.sample,self.unit = wildcards
        self.genome = str(resource_dir/"genome.fasta")
        self.threads = threads
        self.outdir = Path(outdir)

        # input
        self.trimmed_reads1,self.trimmed_reads2 = self.outdir/f"trimmed/{self.sample}-{self.unit}.1.fastq.gz",self.outdir/f"trimmed/{self.sample}-{self.unit}.2.fastq.gz"
        # output
        mapped_dir = self.outdir/"mapped"
        mapped_dir.mkdir(parents=True,exist_ok=True)
        
        self.output_bam = self.outdir/f"mapped/{self.sample}-{self.unit}.sorted.bam"
        # log
        self.log_dir = self.outdir/"logs/mapped"
        self.log_dir.mkdir(parents=True,exist_ok=True)
        

    def map_reads(self):

        _logs = self.log_dir/f"{self.sample}-{self.unit}_bwa.log"
        extra = get_read_group(units,self.sample,self.unit)

        # Extract arguments.
        sort = 'samtools'
        sort_order =  "coordinate"
        sort_extra = ''
        #sort_extra = snakemake.params.get("sort_extra", "")


        # Check inputs/arguments.
        if not isinstance(self.trimmed_reads1, str) and isinstance(self.trimmed_reads2, str) and len([self.trimmed_reads1,self.trimmed_reads2]) not in {
            1,
            2,
        }:
            raise ValueError("input must have 1 (single-end) or " "2 (paired-end) elements")

        if sort_order not in {"coordinate", "queryname"}:
            raise ValueError("Unexpected value for sort_order ({})".format(sort_order))

        # Determine which pipe command to use for converting to bam or sorting.
        if sort == "none":
            # Simply convert to bam using samtools view.
            pipe_cmd = f"samtools view -Sbh -o {self.output_bam} -"
        elif sort == "samtools":
            # Sort alignments using samtools sort.
            pipe_cmd = f"samtools sort {sort_extra} -o {self.output_bam} -"
            # Add name flag if needed.
            if sort_order == "queryname":
                sort_extra += " -n"
            prefix = os.path.splitext(self.output_bam)[0]
            sort_extra += " -T " + prefix + ".tmp"
        elif sort == "picard":
            # Sort alignments using picard SortSam.
            pipe_cmd = (
                f"picard SortSam {sort_extra} INPUT=/dev/stdin"
                f" OUTPUT={self.output_bam} SORT_ORDER={sort_order}"
            )
        else:
            raise ValueError("Unexpected value for params.sort ({})".format(sort))

        cmd = (
            "(bwa mem"
            f" -t {self.threads}"
            f" {extra}"
            f" {self.genome}"
            f" {self.trimmed_reads1} {self.trimmed_reads2}"
            " | " + pipe_cmd + f")"
        )
        if not Path(self.outdir/f"mapped/{self.sample}-{self.unit}.sorted.bam").exists():
            debug_subprocess_call(cmd)
        
        
        


class Get_recal():
    """
    get recal bam
    """
    def __init__(self,wildcards,outdir):
        self.sample,self.unit = wildcards
        self.outdir = Path(outdir)

        # input
        self.mapped_bam = self.outdir/f"mapped/{self.sample}-{self.unit}.sorted.bam"
        self.ref = str(resource_dir/"genome.fasta")
        self.refflat = str(resource_dir/"genome.refFlat")
        self.dict = str(resource_dir/"genome.dict")
        self.known = str(resource_dir/"variation.noiupac.vcf.gz")
        self.known_idx = str(resource_dir/"variation.noiupac.vcf.gz.tbi")
        
        # output
        dedup_dir = self.outdir/"dedup"
        recal_dir = self.outdir/"recal"
        qc_dir = self.outdir/"qc/dedup"
        
        dedup_dir.mkdir(parents=True,exist_ok=True)
        recal_dir.mkdir(parents=True,exist_ok=True)
        qc_dir.mkdir(parents=True,exist_ok=True)
        
        ## markduplicates
        self.markduplicates_bam = self.outdir/f"dedup/{self.sample}-{self.unit}.bam"
        self.markduplicates_metrics = self.outdir/f"qc/dedup/{self.sample}-{self.unit}.metrics.txt"
        ## recalibrate_base_qualities
        self.recal_grp = self.outdir/f"recal/{self.sample}-{self.unit}.grp"
        ## apply_base_quality_recal
        self.recal_bam = self.outdir/f"recal/{self.sample}-{self.unit}.bam"
        self.recal_bam_bai = self.outdir/f"recal/{self.sample}-{self.unit}.bam.bai"

        # log
        self.log_dir = self.outdir/"logs/mapped"
        self.log_dir.mkdir(parents=True,exist_ok=True)


    def markduplicates(self):
        _log = self.log_dir/f"MarkDuplicates-{self.sample}-{self.unit}.log"
        extra = config["params"]["picard"]["MarkDuplicates"]
        cmd = (
                "picard MarkDuplicates -Xmx8g "
                f"{extra} "  # Tool and its subcommand
                f"I={self.mapped_bam} "  # Input bam(s)
                f"O={self.markduplicates_bam} "  # Output bam
                f"M={self.markduplicates_metrics} "
                "ASO=coordinate > {_log}"  # Output metrics
            )
        if not Path(self.markduplicates_bam).exists():
            debug_subprocess_call(cmd)


    def recalibrate_base_qual(self):
        """
        
        """
        _log = self.log_dir/f"BaseRecalibrator-{self.sample}-{self.unit}.log"
        bam = str(self.outdir/get_recal_input(self.sample,self.unit))
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
        
        if self.known:
            if isinstance(self.known, str):
                known = self.known
                known = f"--known-sites {known}"

        cmd = (
            f"gatk --java-options '-Xmx8g' BaseRecalibrator {extra} "
            "--cloud-prefetch-buffer 400 "
            f"-R {self.ref} -I {bam} "
            f"-O {self.recal_grp} {known} > {_log}"
        )
        if not Path(self.recal_bam).exists():
            debug_subprocess_call(cmd)
        

    def apply_base_quality_recal(self):
        """
        """
        _log = self.log_dir/f"ApplyBQSR-{self.sample}-{self.unit}.log"
        bam = str(self.outdir/get_recal_input(self.sample,self.unit))
        extra = get_regions_param()
        
        cmd = (
            f"gatk --java-options '-Xmx8g' ApplyBQSR {extra} "
            f"-R {self.ref} -I {bam} "
            f"--bqsr-recal-file {self.recal_grp} "
            f"-O {self.recal_bam} > {_log}" 
        )
        _cmd = (f"samtools index {self.recal_bam} {self.recal_bam_bai}")
        if not Path(self.recal_bam).exists():
            debug_subprocess_call(cmd)
            debug_subprocess_call(_cmd)


    def qc(self):
        _bam = str(self.mapped_bam)
        name = _bam.rsplit('/',1)[1].rsplit('.',2)[0]
        
        output_samtools=f"{str(self.outdir)}/mapped/{name}_mapping.txt"
        cmd_samtools = (f"samtools stats {_bam} > {output_samtools}")
        cmd_python = (f"python {str(CU_PATH.parents[1])}/mapping/collect_feature.py --bam {_bam} --outdir {str(self.outdir)}/mapped --refflat {self.refflat} --sample {name}")
        
        if not Path(f"{str(self.outdir)}/mapped/{name}_mapping.txt").exists():
            debug_subprocess_call(cmd_samtools)
            debug_subprocess_call(cmd_python)

        # coverage
        cmd_coverage = (f'/SGRNJ01/Public/Software/conda_env/pacbio_v1/bin/samtools coverage {str(self.mapped_bam)} > {str(self.outdir)}/mapped/{name}_frag.cov')
        debug_subprocess_call(cmd_coverage)

        with open(f"{str(self.outdir)}/mapped/{name}_frag.cov") as cov:
            i = 0 
            for line in cov.readlines():
                if line.startswith("202"):
                    i += 1

        df = pd.read_csv(f'{str(self.outdir)}/mapped/{name}_frag.cov',skiprows=i,sep="\t")
        # 
        try:
            index = df.loc[df['#rname'] =="#rname"].index.tolist()[0]
            df = df.loc[:index-1]
        except IndexError:
            df = df
        df['covbases'] = df['covbases'].astype("int64")
        df['endpos'] = df['endpos'].astype("int64")
        mean_coverage = df['covbases'].sum()/df['endpos'].sum()
        df['meandepth'] = df['meandepth'].astype("float32")
        df['bases_num'] = df['meandepth']*df['endpos']
        mean_depth = df['bases_num'].sum()/df['endpos'].sum()

        df = df[['#rname','coverage','meandepth']]
        df['coverage'] = df['coverage'].map(lambda x:str(x)[:5]+"%")
        df['meandepth'] = df['meandepth'].map(lambda x:str(x)+"X")
        df.to_csv(f"{str(self.outdir)}/mapped/{name}_chrom.coverage.txt",sep='\t',index=None)
        with open(f"{str(self.outdir)}/mapped/{name}_chrom.coverage.txt",'r+') as f:
            content = f.read()
            f.seek(0,0)
            f.write(f'# mean coverage: {mean_coverage*100}%\n# mean depth: {mean_depth}X\n'+content)

        cmd_clean = (f'rm {str(self.outdir)}/mapped/{name}_frag.cov')
        debug_subprocess_call(cmd_clean)
        
        
    def get_recal(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.markduplicates()
        self.recalibrate_base_qual()
        self.apply_base_quality_recal()
        self.qc()
        self.merge()
        
        
    @add_log
    def merge(self):
        # merge stat
        bam_list = get_bam_file(f"{str(self.outdir)}/mapped") 
        dic = defaultdict(lambda:defaultdict(list))
        for bam in bam_list:
            name = str(bam).rsplit('/',1)[1].rsplit('.',2)[0]
            # raw reads and mapping read ratio
            with pysam.FastxFile(f'{str(self.outdir)}/trimmed/{name}.1.fastq.gz') as fq:
                raw_reads = 0
                for _ in fq:
                    raw_reads += 1
                raw_reads = raw_reads*2
                dic[name]["Raw Reads"] = raw_reads
            with open(f"{str(self.outdir)}/mapped/{name}_mapping.txt") as mapping:
                i = 0
                for line in mapping.readlines():
                    i += 1
                    if (i == 15):
                        mapped_paired_reads = line.split(":")[1].split("#")[0].strip()
                        dic[name]["Mapped Paired Reads"] = mapped_paired_reads
                        break
            mapped_rate = int(mapped_paired_reads)/int(raw_reads)
            dic[name]["Mapped Reads Rate"] = f'{str(mapped_rate*100)[:5]}%'
            # percent duplication
            with open(f"{str(self.outdir)}/qc/dedup/{name}.metrics.txt") as dedup:
                i = 0
                for line in dedup.readlines():
                    i += 1
                    if (i == 8):
                        Duplication = float(line.rsplit("\t",2)[1])
                        dic[name]['Percent Duplication'] = f'{str(Duplication*100)[:5]}%'
                        break
            # mean depth
            with open(f"{str(self.outdir)}/mapped/{name}_chrom.coverage.txt") as cov:
                    for line in cov.readlines():
                        if line.startswith("# mean coverage"):
                            mean_coverage = line.strip().split(":")[1]
                            dic[name]["Mean Coverage"] = mean_coverage
                        if line.startswith("# mean depth"):
                            mean_depth = line.strip().split(":")[1]
                            dic[name]["Mean Depth"] = mean_depth
            # mapped regions
            with open(f"{str(self.outdir)}/mapped/{name}_stat.txt") as stat:
                line_list = [line.strip() for line in stat.readlines()]
                Base_Pairs_Mapped_to_Exonic_Regions = line_list[0].split(":")[1]
                dic[name]["Base Pairs Mapped to Exonic Regions"] = Base_Pairs_Mapped_to_Exonic_Regions
                Base_Pairs_Mapped_to_Intronic_Regions = line_list[1].split(":")[1]
                dic[name]["Base Pairs Mapped to Intronic Regions"] = Base_Pairs_Mapped_to_Intronic_Regions
                Base_Pairs_Mapped_to_Intergenic_Regions = line_list[2].split(":")[1]
                dic[name]["Base Pairs Mapped to Intergenic Regions"] = Base_Pairs_Mapped_to_Intergenic_Regions 
        merge_df = pd.DataFrame.from_dict(dic,orient="index")
        merge_df.to_csv(f"{str(self.outdir)}/mapped/merge.tsv",sep="\t")


@add_log
def run(params):
    """
    """
    wildcards,threads,outdir = params
    app_map = Map_reads(wildcards,threads,outdir)
    app_map.map_reads()
    app_getrecal = Get_recal(wildcards,outdir)
    app_getrecal.get_recal()



def main():
    outdir = config['outdir']
    threads = 8

    param_list = []
    for wildcards in units.index: 
        param_list.append((wildcards,threads,outdir))

    with multiprocessing.Pool(len(param_list)) as p:
        r=list(tqdm(p.map(run,param_list),total=len(param_list),desc='Mapping reads '))
    p.close()
    p.join()

    # clean
    #clean_cmd = (f"rm -rf {outdir}/mapped {outdir}/dedup")
    #debug_subprocess_call(clean_cmd)

if __name__ == '__main__':
    main()