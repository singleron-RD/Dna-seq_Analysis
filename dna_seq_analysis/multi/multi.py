import os
from collections import defaultdict
from pathlib import Path
import importlib


from dna_seq_analysis.__init__ import __STEPS__ ,RESOURCE
from dna_seq_analysis.tools.common import *


CP = Path(__file__).absolute()
ROOT_PATH = str(CP.parents[1])

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass



def find_step_module(step):
    file_path_dict = {
        'step': f'{ROOT_PATH}/{step}/{step}.py'
    }

    if os.path.exists(file_path_dict['step']):
        step_module = importlib.import_module(f"dna_seq_analysis.{step}.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {step}.{step}")

    return step_module



class Multi():
    def __init__(self):
        self.__STEPS__ = __STEPS__
        self.__CONDA__ = os.path.basename(os.environ['CONDA_DEFAULT_ENV'])
        self.__APP__ = 'dna-seq'
        self.steps_not_run = ['ref']

        # remove
        for step in self.steps_not_run:
            if step in self.__STEPS__:
                self.__STEPS__.remove(step)

        # add args
        self.parser = None
        self.common_args()
        self.step_args()

        # set
        self.args = None
        self.col4_default = None
        self.last_step = ''
        self.fq_suffix = ""
        self.steps_run = self.__STEPS__
        self.logdir = None
        self.sjm_dir = None
        self.sjm_file = None

        self.sjm_cmd = ''
        self.sjm_order = ''
        self.shell_dict = defaultdict(str)

        self.outdir_dic = {}

    def common_args(self):
        readme = f'{self.__APP__} multi-samples'
        parser = argparse.ArgumentParser(readme,
                                         formatter_class=ArgFormatter,
                                         conflict_handler='resolve')
        parser.add_argument('--mod', help='Which type of script to generate, `sjm` or `shell`.',
            choices=['sjm', 'shell'], default='sjm')
        parser.add_argument('--queue', help='Only works if the `--mod` selects `sjm`.')
        parser.add_argument('--rm_files', action='store_true',
            help='Remove redundant fastq,bam,vcf files after running.')
        parser.add_argument('--steps_run', 
            help='''
Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `split` and `trimming`, 
use `--steps_run split,trimming`
''', 
            default='all')
        # sub_program parser do not have
        parser.add_argument('--outdir', help='Output directory.', default="./")
        parser.add_argument('--thread', help='Number of threads', default=4)
        parser.add_argument('--memory', help='The memory size used, the default is GB', default=10)
        self.parser = parser
        return parser


    def step_args(self):
        for step in self.__STEPS__:
            step_module = find_step_module(step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(self.parser, sub_program=False)


    def prepare(self):
        """
        parse_mapfile, make log dir, init script variables, init outdir_dic
        make sjm dir, sjm file
        """
        self.args = self.parser.parse_args()

        if self.args.steps_run != 'all':
            self.steps_run = self.args.steps_run.strip().split(',')
        
        if self.args.mod == 'sjm':

            self.sjm_dir = Path(f'{self.outdir}/sjm/')
            self.sjm_dir.mkdir(parents=True,exist_ok=True)
            self.logdir = Path(self.args.outdir + '/log')
            self.logdir.mkdir(parents=True,exist_ok=True)

            self.sjm_file = f'{self.sjm_dir}/sjm.job'
            self.sjm_cmd = f'log_dir {self.logdir}\n'

        index = 0
        for step in self.__STEPS__:
            step_outdir = f"{self.args.outdir}/{index:02d}.{step}"
            self.outdir_dic.update({step: step_outdir})
            index += 1

    def generate_cmd(self, cmd, step, sample, m=1, x=1):
        if sample:
            sample = "_" + sample
        sched_options = f'sched_options -w n -cwd -V -l vf={m}g,p={x}'
        if self.args.queue:
            sched_options += f' -q {self.args.queue} '
        self.sjm_cmd += f'''
job_begin
    name {step}{sample}
    {sched_options}
    cmd source activate {self.__CONDA__}; {cmd}
job_end
'''

    def process_cmd(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=m, x=x)
        self.shell_dict[sample] += cmd + '\n'
        if self.last_step:
            self.sjm_order += f'order {step}_{sample} after {self.last_step}_{sample}\n'
        self.last_step = step


    def parse_step_args(self, step):
        step_module = find_step_module(step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        step_parser = argparse.ArgumentParser(step_module)
        func_opts(step_parser, sub_program=False)
        args = step_parser.parse_known_args()
        return args


    def get_cmd_line(self, step):
        """ get cmd line without input
            return str
        """
        args = self.parse_step_args(step)
        args_dict = args[0].__dict__
        step_prefix = (
            f'{self.__APP__}_{step} '
            f'--config_path {self.outdir_dic[step]} '
            f'--outdir {self.outdir_dic[step]} '

        )
        cmd_line = step_prefix

        for arg in args_dict:
            if args_dict[arg] is False:
                continue
            if args_dict[arg] is True:
                cmd_line += f'--{arg} '
            else:
                if args_dict[arg]:
                    matches = [' ', '-']
                    arg_string = str(args_dict[arg])
                    if any(char in arg_string for char in matches):  # need quote
                        cmd_line += f'--{arg} "{arg_string}" '
                    else:
                        cmd_line += f'--{arg} {arg_string} '


    def ref(self):
        step = "ref" 
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '   
        )
        self.process_cmd(cmd, step, m=1, x=1)

    

    def run_steps(self):
        for sample in self.fq_dict:
            self.last_step = ''
            for step in self.steps_run:
                if step in self.steps_not_run:
                    continue
                try:
                    method_to_call = getattr(self, step)
                except AttributeError as attr_not_exist:
                    raise NotImplementedError(
                        "Class `{}` does not implement `{}`".format(self.__class__.__name__, step)
                    ) from attr_not_exist
                method_to_call(sample)


    def end(self):
        if self.args.mod == 'sjm':
            with open(self.sjm_file, 'w') as fh:
                fh.write(self.sjm_cmd + '\n')
                fh.write(self.sjm_order)
        if self.args.mod == 'shell':
            os.system('mkdir -p ./shell/')
            for sample in self.shell_dict:
                with open(f'./shell/{sample}.sh', 'w') as f:
                    f.write(self.shell_dict[sample])

    def run(self):
        self.prepare()
        self.run_steps()
        self.end()



def main():
    outdir = config['outdir']
    pass