import pandas as pd
import methlab as ml
import pytest

mate1 = [
    '/long/path/282462_TTCTCGTGCAATACACAGAG_S23_R1_001.fastq.gz',
    '/long/path/282462_TCGTGCATTCGCGTTGGTAT_S5_R1_001.fastq.gz',
    '/long/path/282462_GCCTAACGTGAGCTCTCAAG_S62_R1_001.fastq.gz',
    '/long/path/282462_CATTCACGCTGTTGTACTCA_S26_R1_001.fastq.gz',
    '/long/path/282462_GCCATATAACACACAATATC_S81_R1_001.fastq.gz',
    '/long/path/282462_AACCAGCCACGCCACAGCAC_S60_R1_001.fastq.gz'
]
mate2 = [
    '/long/path/282462_TTCTCGTGCAATACACAGAG_S23_R2_001.fastq.gz',
    '/long/path/282462_TCGTGCATTCGCGTTGGTAT_S5_R2_001.fastq.gz',
    '/long/path/282462_GCCTAACGTGAGCTCTCAAG_S62_R2_001.fastq.gz',
    '/long/path/282462_CATTCACGCTGTTGTACTCA_S26_R2_001.fastq.gz',
    '/long/path/282462_GCCATATAACACACAATATC_S81_R2_001.fastq.gz',
    '/long/path/282462_AACCAGCCACGCCACAGCAC_S60_R2_001.fastq.gz'
]

# Example of a sample sheet giving indices 
sample_sheet = pd.read_csv("tests/test_data/NGS_sample_sheet.csv")

class Test_align_fastq_with_plate_positions:

    def test_indices_match_filenames(self):
        """
        Test that align_fastq_with_plate_positions gives the correct fastq files.
        """
        new_sheet = ml.align_fastq_with_plate_positions(mate1, mate2, sample_sheet)
        # Check the indices appear in the filenames
        assert all(
            [ new_sheet.iloc[i]['index1'] in new_sheet.iloc[i]['fastq1'] for i in range(new_sheet.shape[0]) ]
        )
        assert all(
            [ new_sheet.iloc[i]['index2'] in new_sheet.iloc[i]['fastq1'] for i in range(new_sheet.shape[0]) ]
        )
    
    def test_correct_order(self):
        """
        confirm that the output gives files in the right order, even if the input is not
        """
        shuffled_mate1 = [mate1[i] for i in [3,4,5,0,1,2]]
        shuffled_mate2 = [mate2[i] for i in [2,5,3,4,1,0]]
        new_sheet = ml.align_fastq_with_plate_positions(shuffled_mate1, shuffled_mate2, sample_sheet)
        # Check the indices appear in the filenames
        assert all(
            [ new_sheet.iloc[i]['index1'] in new_sheet.iloc[i]['fastq1'] for i in range(new_sheet.shape[0]) ]
        )
        assert all(
            [ new_sheet.iloc[i]['index2'] in new_sheet.iloc[i]['fastq1'] for i in range(new_sheet.shape[0]) ]
        )

    def test_lists_are_different_lengths(self):
        with pytest.raises(Exception):
            ml.align_fastq_with_plate_positions(mate1, mate2[:3], sample_sheet)

    def test_lists_empty(self):
        with pytest.raises(Exception):
            ml.align_fastq_with_plate_positions([], [], sample_sheet)

    def test_catch_duplicate_indices(self):
        """
        Test that a warning is raised if one or more pairs of input files contain
        the same indices.
        """
        with pytest.raises(Exception):
            ml.align_fastq_with_plate_positions(
                mate1 + [mate1[0]],
                mate2 + [mate2[0]],
                sample_sheet
                )