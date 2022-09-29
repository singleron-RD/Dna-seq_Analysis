import argparse
import os




class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

class Multi():

    def __init__(self, assay):
        self.__ASSAY__ = assay
        #self.__STEPS__ = init_module.__STEPS__
        self.__CONDA__ = os.path.basename(os.environ['CONDA_DEFAULT_ENV'])
        self.__APP__ = 'dna-seq-analysis'
       

        # add args
        self.parser = None
        self.common_args()
        self.step_args()

        # set
        self.args = None     

    def common_args(self):
        readme = f'{self.__ASSAY__} multi-samples'
        parser = argparse.ArgumentParser(readme,
                                         formatter_class=ArgFormatter,
                                         conflict_handler='resolve')
        parser.add_argument(
            '--config_path',
            help='The position where the config.yaml is located',
            required=True)
        
        #parser.add_argument('--steps_run', 
        #    help='''
#Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `split` and `trimming`, 
#use `--steps_run split,trimming`
#''', 
#            default='all')
        self.parser = parser
        return parser

