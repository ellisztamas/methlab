import pytest
import numpy as np

from methlab.BismarkSam import *

sam = read_SAM("tests/test_data/chloroplast.sam") # Example with no YS tags
strandID=read_SAM("tests/test_data/strandID.sam") # example with YS tags.

def test_read_SAM():
    """
    Check import of a sample sam file.
    """
    # Check they all import properly
    assert all( [isinstance(read, BismarkSam) for read in sam] )
    # Check individual attributes work correctly
    assert isinstance(sam[0].id, str)
    assert isinstance(sam[0].length, int)
    assert sam[0].chr == 'ChrC'
    assert( len(sam[0].seq) == len(sam[0].xm_tag) )
    
def test_counting_functions():
    assert sam[0].count_mC() == [0,32]
    # Read 0 should have zero methylated Cs, 32 unmethylated, be 135bp and show no clusters
    assert sam[0].mC_per_read() == ["ChrC", 0, 32, len(sam[0].xm_tag), 'NA']
    # Read 2 should have 15 methylated Cs, 16 unmethylated, be 135bp and does show clusters
    assert sam[2].mC_per_read() == ['ChrC', 15,6, len(sam[2].xm_tag), True]

def test_strand():
    # This SAM file should have only NA entries for strand because Bismark was run without the --strandID flag
    assert all([read.strand == 'NA' for read in sam])
    # Entries in strandID should all be one of 'OT', "OB", 'CTOT', 'CTOB'
    assert strandID[0].strand == 'CTOB'
    all([read.strand in ['OT', "OB", 'CTOT', 'CTOB'] for read in strandID])
