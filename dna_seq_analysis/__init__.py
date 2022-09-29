__VERSION__ = "1.0"
__version__ = __VERSION__
__all__ = ['tools']

__STEPS__ = [
    'ref',
    'split',
    'trimming',
    'mapping',
    'calling',
    'filtering',
    'annotation',
    'stats',
    'qc',
    ]

__ASSAY__ = 'dna-seq'

RESOURCE = {
    'ref': {'m': 1, 'x': 1},
    'split': {'m': 5, 'x': 1},
    'trimming': {'m': 5, 'x': 1},
    'mapping': {'m': 30, 'x': 4},
    'calling': {'m': 20, 'x': 8},
    'filtering': {'m': 5, 'x': 1},
    'annotation': {'m': 5, 'x': 1},
    'stats': {'m': 1, 'x': 1},
    'qc': {'m': 5, 'x': 1},
}