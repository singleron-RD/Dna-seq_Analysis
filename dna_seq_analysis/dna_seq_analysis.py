import argparse
import importlib
import os

from dna_seq_analysis.__init__ import __VERSION__,__STEPS__,ROOT_PATH


def find_step_init(step):
    init_module = importlib.import_module(f"dna_seq_analysis.{step}.__init__")
    return init_module


def find_step_module(step):
    file_path_dict = {
        'step': f'{ROOT_PATH}/{step}/{step}.py',
    }
    if os.path.exists(file_path_dict['step']):
        step_module = importlib.import_module(f"dna_seq_analysis.{step}.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {step}.{step}")

    return step_module



class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main():
    """dna seq analysis cli
    """
    parser = argparse.ArgumentParser(description='WGS process', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers(dest='subparser_assay')
    
    for step in __STEPS__:
        # import function and opts
        step_module = find_step_module(step)
        func = getattr(step_module, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        parser_step = subparsers.add_parser(step, description=f'{step.upper()} of WGS process',formatter_class=ArgFormatter)
        func_opts(parser_step)
        parser_step.set_defaults(func=func)

    args = parser.parse_args()
    if len(args.__dict__) <= 1:
        # No arguments or subcommands were given.
        parser.print_help()
        parser.exit()
    else:
        args.func(args)

    return args

if __name__ == '__main__':
    main()