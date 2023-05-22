import os
from pathlib import Path
import argparse

from dna_seq_analysis.tools.common import *
from dna_seq_analysis.__init__ import __STEPS__

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


class Multi():
    def __init__(self):

        self.__STEPS__ = __STEPS__
        self.__CONDA__ = os.path.basename(os.environ['CONDA_DEFAULT_ENV'])
        self.__APP__ = 'dna-seq'
    
        # set
        self.last_step = ''
        self.steps_run = self.__STEPS__
        self.logdir = None
        self.sjm_dir = None
        self.sjm_file = None

        self.sjm_cmd = ''
        self.sjm_order = ''
        self.shell = ''
        
        # add args
        self.parser = None
        self.common_args()
        self.step_args()
        
        self.steps_not_run = []
        
    def common_args(self):
        readme = f'WGS process multi-samples'
        parser = argparse.ArgumentParser(readme,
                                         formatter_class=ArgFormatter,
                                         conflict_handler='resolve')
        parser.add_argument('--mod', help='Which type of script to generate, `sjm` or `shell`.',
            choices=['sjm', 'shell'], default='sjm')
        parser.add_argument('--queue', help='Only works if the `--mod` selects `sjm`.')
        parser.add_argument('--rm_files', action='store_true',
            help='Remove redundant fastq and bam files after running.')
        parser.add_argument('--steps_run', 
            help='''
    Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `split` and `trimming`, 
    use `--steps_run split,trimming`
    ''', 
            default='all')
        parser.add_argument('--steps_not_run', 
            help='''
    Steps not to run. Multiple Steps are separated by comma. For example, if you only want not to run `ref` and `split`, 
    use `--steps_not_run ref,split`
    ''')
        parser.add_argument('--whether_split', help='Whether to split the original sequence according to the index.',choices=['true','false'],required=True)
        # sub_program parser do not have
        parser.add_argument('--thread', help='Number of threads.', default=4,type=int)
        parser.add_argument('--outdir', help='The directory of script output.', default='./')
        parser.add_argument('--memory', help='The memory size used, the default is GB', default=10,type=int)
        self.parser = parser
        return parser


    def step_args(self):
        for step in self.__STEPS__:
            step_module = find_step_module(step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(self.parser, sub_program=True)
    

    def prepare(self):
        """
        parse_config, make log dir, init script variables, init outdir_dic
        make sjm dir, sjm file
        """
        self.args = self.parser.parse_args()

        if self.args.steps_run != 'all':
            self.steps_run = self.args.steps_run.strip().split(',')
            self.steps_not_run = [step for step in self.__STEPS__ if step not in self.steps_run]
        
        if self.args.whether_split == 'false':
            self.steps_not_run += ['split']
            
        if self.args.whether_split == 'true':
            if not self.args.mapfile:
                raise ValueError(f"If need to split the original sequence according to the index, please use `--mapfile` parameter provide the corresponding mapfile file.")
        
        if self.args.steps_not_run :
            self.steps_not_run += self.args.steps_not_run.strip().split(',')
            # remove
            for step in self.steps_not_run:
                if step in self.__STEPS__:
                    self.__STEPS__.remove(step)
        
        if 'ref' not in self.steps_not_run:
            if not (self.args.species and self.args.release and self.args.build):
                raise ValueError(f"Please use `--species --release --build` parameter provide correct genome information.")
            
        if 'annotation' not in self.steps_not_run:
            if not (self.args.species and self.args.build):
                raise ValueError(f"Please use `--species --build` parameter provide correct ensembl species name and genome build.")
        
        if self.args.mod == 'sjm':

            self.sjm_dir = f'{self.args.outdir}/sjm/'
            Path(self.sjm_dir).mkdir(parents=True,exist_ok=True)
            self.logdir = self.args.outdir + '/log'
            Path(self.logdir).mkdir(parents=True,exist_ok=True)

            self.sjm_file = f'{self.sjm_dir}/sjm.job'
            self.sjm_cmd = f'log_dir {self.logdir}\n'

        # parse_config
        self.config = parse_config(self.args.config_path)
        
        # thread and memory
        self.thread = self.args.thread
        self.memory = self.args.memory

    
    def generate_cmd(self, cmd, step, m=1, x=1):
    
        sched_options = f'sched_options -w n -cwd -V -l vf={m}g,p={x}'
        if self.args.queue:
            sched_options += f' -q {self.args.queue} '
        self.sjm_cmd += f'''
job_begin
    name {step}
    {sched_options}
    cmd source activate {self.__CONDA__}; {cmd}
job_end
'''

    def process_cmd(self, cmd, step,  m=1, x=1):
        self.generate_cmd(cmd, step,  m=m, x=x)
        self.shell += cmd + '\n'
        if self.last_step:
            self.sjm_order += f'order {step} after {self.last_step}\n'
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
            f'{self.__APP__} {step} '
            f'--config_path {self.args.config_path} '
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

        return cmd_line


    def ref(self):
        step = "ref"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=self.memory, x=1)
        
    def split(self):
        step = "split"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=self.memory, x=1)
        
    def trim(self):
        step = "trim"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=self.memory, x=self.thread)
        
    def map(self):
        step = "map"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=30, x=self.thread)
        
    def call(self):
        step = "call"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=30, x=self.thread)

    def filter(self):
        step = "filter"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=self.memory, x=self.thread)
    
    def annotation(self):
        step = "annotation"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=self.memory, x=self.thread)
    
    def stat(self):
        step = "stat"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=5, x=1)
        
    def qc(self):
        step = "qc"
        cmd_line = self.get_cmd_line(step)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, m=5, x=1)
        

    def run_steps(self):
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
            method_to_call()
                

    def end(self):
        if self.args.mod == 'sjm':
            with open(self.sjm_file, 'w') as fh:
                fh.write(self.sjm_cmd + '\n')
                fh.write(self.sjm_order)
        if self.args.mod == 'shell':
            text_head = (f'source activate {self.__CONDA__}\n'
                #"#!/usr/bin/env bash\n"
                #"set -euxo pipefail\n"
                )
            os.system(f'mkdir -p {self.args.outdir}/shell/')
            with open(f'{self.args.outdir}/shell/run.sh', 'w') as f:
                f.write(text_head)
                f.write(self.shell)

    def run(self):
        self.prepare()
        self.run_steps()
        self.end()


def main():
    runner = Multi()
    runner.run()