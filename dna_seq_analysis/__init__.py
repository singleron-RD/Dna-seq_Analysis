import os

__VERSION__ = "1.0"
__version__ = __VERSION__


__STEPS__ = [
    'ref',
    'split',
    'trim',
    'map',
    'call',
    'filter',
    'annotation',
    'stat',
    'qc'
    ]
__ASSAY__ = 'dna-seq'

ROOT_PATH = os.path.dirname(__file__)