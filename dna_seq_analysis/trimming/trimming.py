from email.policy import default
from pathlib import Path
import multiprocessing
from tqdm import tqdm

from dna_seq_analysis.tools.common import *



# Trim reads
class Trim_reads():
    """
    trim reads
    """
    def __init__(self,args,fastq,wildcards,trimmer,mod):
        self.fastq = fastq
        self.mod = mod
        self.outdir = Path(args.outdir) #
        self.sample,self.unit = wildcards
        self.trimmer = trimmer
        self.threads = int(args.thread) #

        # log
        self.log_dir = self.outdir/"logs/trim_reads"
        self.log_dir.mkdir(parents=True,exist_ok=True)
    
    @staticmethod
    def distribute_threads(input_files, output_files, available_threads):
        """
        Distribute available threads between trimmomatic itself and any potential pigz instances
        """
        gzipped_input_files = sum(1 for file in input_files if file.endswith(".gz"))
        gzipped_output_files = sum(1 for file in output_files if file.endswith(".gz"))
        potential_threads_per_process = available_threads // (
            1 + gzipped_input_files + gzipped_output_files
        )
        if potential_threads_per_process > 0:
            # decompressing pigz creates at most 4 threads
            pigz_input_threads = (
                min(4, potential_threads_per_process) if gzipped_input_files != 0 else 0
            )
            pigz_output_threads = (
                (available_threads - pigz_input_threads * gzipped_input_files)
                // (1 + gzipped_output_files)
                if gzipped_output_files != 0
                else 0
            )
            trimmomatic_threads = (
                available_threads
                - pigz_input_threads * gzipped_input_files
                - pigz_output_threads * gzipped_output_files
            )
        else:
            # not enough threads for pigz
            pigz_input_threads = 0
            pigz_output_threads = 0
            trimmomatic_threads = available_threads
        return trimmomatic_threads, pigz_input_threads, pigz_output_threads

    @staticmethod
    def compose_input_gz(filename, threads):
        if filename.endswith(".gz") and threads > 0:
            return "<(pigz -p {threads} --decompress --stdout {filename})".format(
                threads=threads, filename=filename
            )
        return filename

    @staticmethod
    def compose_output_gz(filename, threads, compression_level):
        if filename.endswith(".gz") and threads > 0:
            return ">(pigz -p {threads} {compression_level} > {filename})".format(
                threads=threads, compression_level=compression_level, filename=filename
            )
        return filename


    def trim_read(self):
        _logs = self.log_dir/f'{self.sample}-{self.unit}.log'
        compression_level = "-5"
        # Distribute threads
        input_files = [self.fastq['r1'], self.fastq['r2']]
        trimmed_dir = self.outdir/"trimmed"
        trimmed_dir.mkdir(parents=True,exist_ok=True)
        output_files = [
            f'{trimmed_dir}/{self.sample}-{self.unit}.1.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.1.unpaired.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.2.fastq.gz',
            f'{trimmed_dir}/{self.sample}-{self.unit}.2.unpaired.fastq.gz',
        ]

        trimmomatic_threads, input_threads, output_threads = self.distribute_threads(
            input_files, output_files, self.threads
        )

        input_r1, input_r2 = [
            self.compose_input_gz(filename, input_threads) for filename in input_files
        ]

        output_r1, output_r1_unp, output_r2, output_r2_unp = [
            self.compose_output_gz(filename, output_threads, compression_level)
            for filename in output_files
        ]
        
        cmd = (
                f"trimmomatic {self.mod} -threads {trimmomatic_threads} "
                f"-trimlog {trimmed_dir}/{self.sample}-{self.unit}.trimlog.txt "
                f"{input_r1} {input_r2} "
                f"{output_r1} {output_r1_unp} "
                f"{output_r2} {output_r2_unp} "
                f"{self.trimmer} > {_logs}"
            )
        debug_subprocess_call(cmd)


@add_log
def run(params):
    """
    """
    args,fastq,wildcards,trimmer,mod = params
    app = Trim_reads(args=args,fastq=fastq,wildcards=wildcards,trimmer=trimmer,mod=mod)
    app.trim_read()


def main():
    #outdir = config['outdir']
    param_list = []
    parser = s_common()
    parser = get_opts_trimming(parser,sub_program=False)
    args = parser.parse_args()

    config_path = args.config_path
    configfile = Path(config_path)/"config.yaml"
    with open(configfile,"r") as test_file:
        config = yaml.load(test_file,Loader=yaml.FullLoader)
    units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)


    for wildcards in units.index:
        fastq = get_fastq(units,wildcards)
        mod = 'PE' if len(fastq) > 1 else 'SE'
        trimmer = " ".join(config["params"]["trimmomatic"][f"{mod}".swapcase()]['trimmer'])
        param_list.append((args,fastq,wildcards,trimmer,mod))

    with multiprocessing.Pool(len(param_list)) as p:
        r=list(tqdm(p.map(run,param_list),total=len(param_list),desc='Trim reads '))
    p.close()
    p.join()


def get_opts_trimming(parser, sub_program=True):
    parser.add_argument(
        '--thread',
        help='Thread to use.',
        default = 4
    )
    if sub_program:
        print('')
    return parser