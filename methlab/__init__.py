"""Top-level package for methlab."""

__author__ = """Tom Ellis"""
__email__ = 'thomas.ellis@gmi.oeaw.ac.at'
__version__ = '0.9.2'

# import argparse

from methlab.align_fastq_with_plate_positions import align_fastq_with_plate_positions
from methlab.CytosineCoverageFile import CytosineCoverageFile
from methlab.BismarkSam import *
from methlab.methylation_state import methylation_state
from methlab.estimate_beta_parameters import estimate_beta_parameters
from methlab.ibdpainting import *