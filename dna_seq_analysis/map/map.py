import os
import multiprocessing
import numpy as np 
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import defaultdict
from pathlib import Path
import unittest
import sys

from dna_seq_analysis.tools.common import *


class Map_reads():
    """
    bwa mapping
    """
    def __init__(self,wildcards,threads,outdir,resource_dir,units,args):
        self.sample,self.unit = wildcards
        self.genome = f'{resource_dir}/genome.fasta'
        self.threads = threads
        self.outdir = Path(outdir)
        self.units = units
        self.args = args
        
        # input
        if self.args.downsample:
            self.trimmed_reads1,self.trimmed_reads2 = self.outdir/f"01.trimmed_ds/{self.sample}-{self.unit}.1.fastq",self.outdir/f"01.trimmed_ds/{self.sample}-{self.unit}.2.fastq"
        else:
            self.trimmed_reads1,self.trimmed_reads2 = self.outdir/f"01.trimmed/{self.sample}-{self.unit}.1.fastq.gz",self.outdir/f"01.trimmed/{self.sample}-{self.unit}.2.fastq.gz"
        # output
        mapped_dir = self.outdir/"02.mapped"
        mapped_dir.mkdir(parents=True,exist_ok=True)
        
        self.output_bam = self.outdir/f"02.mapped/{self.sample}-{self.unit}.sorted.bam"
        # log
        self.log_dir = self.outdir/"logs/mapped"
        self.log_dir.mkdir(parents=True,exist_ok=True)
    
    
    @add_log
    def map_reads(self):
        extra = get_read_group(self.units,self.sample,self.unit)
        # Extract arguments.
        sort = 'samtools'
        sort_order =  "coordinate"
        sort_extra = ''
    
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
            " | " + pipe_cmd + ")"
        )
        cmd_index=(
            f'samtools index {self.output_bam}'
        )
        if not Path(self.outdir/f"02.mapped/{self.sample}-{self.unit}.sorted.bam").exists():
            debug_subprocess_call(cmd)
            debug_subprocess_call(cmd_index)


class Get_recal():
    """
    get recal bam
    """
    def __init__(self,wildcards,threads,outdir,resource_dir,args):
        self.sample,self.unit = wildcards
        self.threads = threads
        self.outdir = Path(outdir)
        self.args = args

        # input
        self.mapped_bam = self.outdir/f"02.mapped/{self.sample}-{self.unit}.sorted.bam"
        self.ref = f"{resource_dir}/genome.fasta"
        self.refflat = f"{resource_dir}/genome.refFlat"
        self.dict = f"{resource_dir}/genome.dict"
        self.known = f"{resource_dir}/variation.noiupac.vcf.gz"
        self.known_idx = f"{resource_dir}/variation.noiupac.vcf.gz.tbi"
        
        # output
        dedup_dir = self.outdir/"03.dedup"
        recal_dir = self.outdir/"04.recal"
        qc_dir = self.outdir/"qc/dedup"
        
        dedup_dir.mkdir(parents=True,exist_ok=True)
        recal_dir.mkdir(parents=True,exist_ok=True)
        qc_dir.mkdir(parents=True,exist_ok=True)
        
        ## markduplicates
        self.markduplicates_bam = self.outdir/f"03.dedup/{self.sample}-{self.unit}.bam"
        self.markduplicates_metrics = self.outdir/f"qc/dedup/{self.sample}-{self.unit}.metrics.txt"
        ## recalibrate_base_qualities
        self.recal_grp = self.outdir/f"04.recal/{self.sample}-{self.unit}.grp"
        ## apply_base_quality_recal
        self.recal_bam = self.outdir/f"04.recal/{self.sample}-{self.unit}.bam"
        self.recal_sorted_bam = self.outdir/f"04.recal/{self.sample}-{self.unit}.sorted.bam"
        self.recal_bam_bai = self.outdir/f"04.recal/{self.sample}-{self.unit}.bam.bai"

        # log
        self.log_dir = self.outdir/"logs/mapped"
        self.log_dir.mkdir(parents=True,exist_ok=True)


    def markduplicates(self):
        _log = self.log_dir/f"MarkDuplicates-{self.sample}-{self.unit}.log"
        if self.args.remove_duplicates:
            extra = 'REMOVE_DUPLICATES=true'
        else:
            extra = ''
        cmd = (
                "picard MarkDuplicates -Xmx8g "
                f"{extra} "  # Tool and its subcommand
                f"I={self.mapped_bam} "  # Input bam(s)
                f"O={self.markduplicates_bam} "  # Output bam
                f"M={self.markduplicates_metrics} > {_log}"  # Output metrics
            )
        cmd_index = (
            f'samtools index {self.markduplicates_bam}'
        )
        if not Path(self.markduplicates_bam).exists():
            debug_subprocess_call(cmd)
            debug_subprocess_call(cmd_index)


    def recalibrate_base_qual(self):
        _log = self.log_dir/f"BaseRecalibrator-{self.sample}-{self.unit}.log"
        
        bam = str(self.outdir/get_recal_input(self.sample,self.unit,self.args.remove_duplicates))
        extra = get_regions_param(self.args.intervals,self.args.interval_padding)
        
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
        bam = str(self.outdir/get_recal_input(self.sample,self.unit,self.args.remove_duplicates))
        extra = get_regions_param(self.args.intervals,self.args.interval_padding)
        
        cmd = (
            f"gatk --java-options '-Xmx8g' ApplyBQSR {extra} "
            f"-R {self.ref} -I {bam} "
            f"--bqsr-recal-file {self.recal_grp} "
            f"-O {self.recal_bam} > {_log}" 
        )
        cmd_sort = (f"samtools sort -@{self.threads} {self.recal_bam} -o {self.recal_sorted_bam}")
        cmd_index = (f"samtools index {self.recal_bam} {self.recal_bam_bai}")
        if not Path(self.recal_bam).exists():
            debug_subprocess_call(cmd)
            debug_subprocess_call(cmd_sort)
            debug_subprocess_call(cmd_index)


    def qc(self):
        _bam = str(self.mapped_bam)
        name = _bam.rsplit('/',1)[1].rsplit('.',2)[0]
        
        output_samtools=f"{str(self.outdir)}/02.mapped/{name}_mapping.txt"
        cmd_samtools = (f"samtools stats {_bam} > {output_samtools}")
        cmd_python = (f"python {str(CU_PATH.parents[1])}/tools/collect_feature.py --bam {_bam} --outdir {str(self.outdir)}/02.mapped --refflat {self.refflat} --sample {name}")
        
        if not Path(f"{str(self.outdir)}/02.mapped/{name}_mapping.txt").exists():
            debug_subprocess_call(cmd_samtools)
        debug_subprocess_call(cmd_python)

        # coverage
        cmd_coverage = (f'/SGRNJ01/Public/Software/conda_env/pacbio_v1/bin/samtools coverage {str(self.mapped_bam)} > {str(self.outdir)}/02.mapped/{name}_frag.cov')
        
        debug_subprocess_call(cmd_coverage)

        with open(f"{str(self.outdir)}/02.mapped/{name}_frag.cov") as cov:
            i = 0 
            for line in cov.readlines():
                if line.startswith("202"):
                    i += 1

        df = pd.read_csv(f'{str(self.outdir)}/02.mapped/{name}_frag.cov',skiprows=i,sep="\t")
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
        df.to_csv(f"{str(self.outdir)}/02.mapped/{name}_chrom.coverage.txt",sep='\t',index=None)
        with open(f"{str(self.outdir)}/02.mapped/{name}_chrom.coverage.txt",'r+') as f:
            content = f.read()
            f.seek(0,0)
            f.write(f'# mean coverage: {mean_coverage*100}%\n# mean depth: {mean_depth}X\n'+content)

        cmd_clean = (f'rm {str(self.outdir)}/02.mapped/{name}_frag.cov')
        debug_subprocess_call(cmd_clean)
    
        
    @add_log  
    def get_recal(self):
        """
        run markduplicates and recalibrate base qual and apply base quality recal
        """
        self.markduplicates()
        self.recalibrate_base_qual()
        self.apply_base_quality_recal()
        self.qc()


class CNV():
    """
    copy number variations
    """
    def __init__(self,threads,outdir,resource_dir):
        self.genome = f"{resource_dir}/genome.fasta"
        self.refflat = f"{resource_dir}/genome.refFlat"
        self.threads = threads
        self.recal_path = Path(f'{outdir}/04.recal')
        
    def cnv(self):
        bed_file = 'access-1kb.bed'
        cmd_bed = (f"cnvkit.py access {self.refflat} -s 100000 -o {bed_file};sed -i '/^[KGMJ]/d' {bed_file}")
        cmd_cnv = ('cnvkit.py batch '
               f'{str(self.recal_path)}/*.sorted.bam -n '
               f'-t {bed_file} --method wgs --fasta {self.genome} --annotate {self.refflat} '
               f'-p {self.threads} --output-dir {str(self.recal_path)}/cnv '
               '--scatter --drop-low-coverage '
        )
        cmd_plot = (f'cnvkit.py heatmap {str(self.recal_path)}/cnv/*.cns -o {str(self.recal_path)}/cnv/Heatmap.pdf')
        debug_subprocess_call(cmd_bed)
        debug_subprocess_call(cmd_cnv)
        debug_subprocess_call(cmd_plot)
        Path(bed_file).unlink()

    def plot_bar(self):
        cnvkit_dir = self.recal_path/'cnv'
        cnr_dict ={}
        for cnr_file in cnvkit_dir.rglob('*sorted.cnr'):
            sample = cnr_file.name.rsplit("-",1)[0]
            cnr_dict.update({sample:str(cnr_file)})
            
        chrom = '12'  # set
        fig,axs = plt.subplots(nrows=len(cnr_dict),ncols=1,figsize=(20,3*len(cnr_dict)))

        for i,sample in enumerate(cnr_dict):
            df = pd.read_csv(cnr_dict[sample],header=0,sep="\t")
            df['chromosome'] = df['chromosome'].astype('str')
            df_sub = (
                        df.query("chromosome == @chrom")
                        .drop(columns=['start','end','gene','depth','weight'])
                    )
            df_sub_pos = df_sub.where(lambda d:d.log2 > 0,0)
            df_sub_pos['chromosome'] = np.arange(df_sub_pos.shape[0])
            df_sub_neg = df_sub.mask(lambda d:d.log2 > 0,0)
            df_sub_neg['chromosome'] = np.arange(df_sub_neg.shape[0])
            lims = df_sub['log2'].min(),df_sub['log2'].max()
        
            axs[i].bar(df_sub_pos['chromosome'],df_sub_pos['log2'],color='red')
            axs[i].set_ylim(lims)
            ax = axs[i].twinx()
            ax.bar(df_sub_neg['chromosome'],df_sub_neg['log2'],color='blue')
            ax.set_ylim(lims)
            ax.set_axis_off()
            #set axs[i]
            for spine in ['top','bottom','left','right']:
                axs[i].spines[spine].set_color('none')
            axs[i].tick_params(bottom=False,top=False,left=False,right=False)
            axs[i].set_xticks([])
            axs[i].set_yticks([])
            # add y lable
            fontdict_lable={'size':10,
                    'color':'k',
                    'family':'serif'}
            axs[i].set_ylabel(f'{sample}',fontdict=fontdict_lable)
        
        # set sub plot
        fig.suptitle('CNVs detected',fontsize=16,fontweight='bold')
        fig.supylabel('Log2(sample/reference)',fontsize=16,fontweight='bold') 
        fig.savefig(f"{cnvkit_dir}/Chorm12_cnr.pdf",dpi=1000,bbox_inches='tight')
        fig.savefig(f"{cnvkit_dir}/Chorm12_cnr.png",dpi=1000,bbox_inches='tight')
    
    
    @add_log
    def run(self):
        self.cnv()
        self.plot_bar() 



@add_log
def run(params):
    """
    """
    wildcards,threads,outdir,resource_dir,units,args = params
    app_map = Map_reads(wildcards,threads,outdir,resource_dir,units,args)
    app_map.map_reads()
    app_getrecal = Get_recal(wildcards,threads,outdir,resource_dir,args)
    app_getrecal.get_recal()


@add_log
def merge(outdir,downsample):
    # merge stat
    bam_list = get_bam_file(f"{outdir}/02.mapped") 
    dic = defaultdict(lambda:defaultdict(list))
    for bam in bam_list:
        name = str(bam).rsplit('/',1)[1].rsplit('.',2)[0]
        # raw reads and mapping read ratio
        fq = f'{outdir}/01.trimmed_ds/{name}.1.fastq' if downsample else f'{outdir}/01.trimmed/{name}.1.fastq.gz'
        _ = get_fq_reads_num(fq)
        raw_reads = _*2
        dic[name]["Raw Reads"] = raw_reads
        with open(f"{outdir}/02.mapped/{name}_mapping.txt") as mapping:
            i = 0
            for line in mapping.readlines():
                i += 1
                if (i == 15):
                    mapped_paired_reads = line.split(":")[1].split("#")[0].strip()
                    dic[name]["Mapped Paired Reads"] = mapped_paired_reads
                    break
        mapped_rate = int(mapped_paired_reads)/int(raw_reads)
        # The value is too small to use scientific notation
        dic[name]["Mapped Reads Rate"] = f'{mapped_rate*100}%'
        # percent duplication
        with open(f"{outdir}/qc/dedup/{name}.metrics.txt") as dedup:
            i = 0
            for line in dedup.readlines():
                i += 1
                if (i == 8):
                    Duplication = float(line.rsplit("\t",2)[1])
                    dic[name]['Percent Duplication'] = f'{Duplication*100}%'
                    break
        # mean depth
        with open(f"{outdir}/02.mapped/{name}_chrom.coverage.txt") as cov:
                for line in cov.readlines():
                    if line.startswith("# mean coverage"):
                        mean_coverage = line.strip().split(":")[1]
                        dic[name]["Mean Coverage"] = mean_coverage
                    if line.startswith("# mean depth"):
                        mean_depth = line.strip().split(":")[1]
                        dic[name]["Mean Depth"] = mean_depth
        # mapped regions
        with open(f"{outdir}/02.mapped/{name}_stat.txt") as stat:
            line_list = [line.strip() for line in stat.readlines()]
            Base_Pairs_Mapped_to_Exonic_Regions = line_list[0].split(":")[1]
            dic[name]["Base Pairs Mapped to Exonic Regions"] = Base_Pairs_Mapped_to_Exonic_Regions
            Base_Pairs_Mapped_to_Intronic_Regions = line_list[1].split(":")[1]
            dic[name]["Base Pairs Mapped to Intronic Regions"] = Base_Pairs_Mapped_to_Intronic_Regions
            Base_Pairs_Mapped_to_Intergenic_Regions = line_list[2].split(":")[1]
            dic[name]["Base Pairs Mapped to Intergenic Regions"] = Base_Pairs_Mapped_to_Intergenic_Regions 
    merge_df = pd.DataFrame.from_dict(dic,orient="index")
    merge_df.to_csv(f"{outdir}/02.mapped/merge.tsv",sep="\t")


def downsample_bam(bam,outdir,sample_name,p):
    """
    downsample to the same level
    """
    print('run downsample bam')
    Path(f'{outdir}/02.mapped_ds').mkdir(parents=True,exist_ok=True)
    cmd_downsample = (f'picard DownsampleSam '
                    f'I={bam} O={outdir}/02.mapped_ds/{sample_name}.bam '
                    f'P={p} R=100 ACCURACY=0.00001 STRATEGY=ConstantMemory'
                    )
    cmd_sort = (f'samtools sort -@8 {outdir}/02.mapped_ds/{sample_name}.bam -o {outdir}/02.mapped_ds/{sample_name}.sorted.bam')
    cmd_index = (f'samtools index {outdir}/02.mapped_ds/{sample_name}.sorted.bam')
    debug_subprocess_call(cmd_downsample)
    debug_subprocess_call(cmd_sort)
    debug_subprocess_call(cmd_index)


def run_downsample_bam(param):
    bam,outdir,sample_name,p = param
    downsample_bam(bam,outdir,sample_name,p)


def downsample_fastq(fastq,outdir,fq_name,min_num):
    """
    downsample to the same level
    """
    print('run downsample fastq')
    soft = 'zcat' if str(fastq).endswith('gz') else 'cat'
    lines_num = min_num*4
    cmd = (f'{soft} {fastq}|head -n {lines_num} > {outdir}/01.trimmed_ds/{fq_name}')
    debug_subprocess_call(cmd)


def run_downsample_fastq(param):
    fastq,outdir,fq_name,min_num = param
    downsample_fastq(fastq,outdir,fq_name,min_num)




@add_log
def map(args):
    config_path = args.config_path
    threads = args.thread
    config = parse_config(config_path)
    outdir = config['outdir']
    resource_dir = config['genomedir']
    
    units_file = f'{config_path}/units.tsv'
    units = get_units(units_file)

    if args.downsample:
        fq_list = []
        for fastq in Path(f'{outdir}/01.trimmed').rglob("*[0-9].fastq.gz"):
            fq_list.append(str(fastq))
        with multiprocessing.Pool() as p:
            fq_num = p.map(get_fq_reads_num,fq_list)
        min_num = min(fq_num)
        Path(f'{outdir}/01.trimmed_ds').mkdir(parents=True,exist_ok=True)
        p_list = []
        for fastq in Path(f'{outdir}/01.trimmed').rglob("*[0-9].fastq.gz"):
            fq_name = fastq.name.rsplit(".",1)[0]
            p_list.append((fastq,outdir,fq_name,min_num))
        with multiprocessing.Pool(len(p_list)) as p:
            p.map(run_downsample_fastq,p_list)
        p.close()
        p.join() 


    param_list = []
    for wildcards in units.index: 
        param_list.append((wildcards,threads,outdir,resource_dir,units,args))

    # Limit the number of CPUs used
    if len(param_list) > 20:
        num_cpu = 20
    else:
        num_cpu = len(param_list)
    
    with multiprocessing.Pool(num_cpu) as p:
        list(tqdm(p.imap(run,param_list),total=len(param_list),unit_scale = True,ncols = 70,file = sys.stdout,desc='Mapping reads '))
    p.close()
    p.join()

    # merge data
    merge(str(outdir),args.downsample)
    
    # CNV
    #cnv_app = CNV(threads,outdir,resource_dir)
    #cnv_app.run()
    
    # clean
    clean_cmd = (f"rm -rf {outdir}/02.mapped {outdir}/03.dedup")
    if args.rm_files:
        debug_subprocess_call(clean_cmd)
    
    

def get_opts_map(parser, sub_program=True):
    parser.add_argument('--thread',help='Number of threads.', default=4,type=int)
    parser.add_argument('--remove_duplicates',help='Delete the duplicate sequence after comparison.', action='store_true')
    parser.add_argument('--intervals',help='One or more genomic intervals over which to operate.This argument may be specified 0 or more times.')
    parser.add_argument('--interval_padding',help='Amount of padding (in bp) to add to each interval you are including.',default=0,type=int)
    parser.add_argument('--rm_files',help='Remove redundant bam files after running.',action='store_true')
    parser.add_argument('--downsample',help='Downsample the raw reads of all samples to a consistent level.',action='store_true')
    
    if sub_program:
        parser = s_common(parser)
    return parser


if __name__ == '__main__':
    unittest.main()