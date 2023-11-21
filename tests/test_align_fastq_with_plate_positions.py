import pandas as pd

import methlab as ml
file_list = [
    'H3H7YDRXY_1#144456_ACTCGCTAAAGGCTAT.fastq.gz',
    'H3H7YDRXY_1#144456_ACTCGCTACCTAGAGT.fastq.gz',
    'H3H7YDRXY_1#144456_ACTCGCTACGTCTAAT.fastq.gz',
    'H3H7YDRXY_2#144456_ACTCGCTAAAGGCTAT.fastq.gz',
    'H3H7YDRXY_2#144456_ACTCGCTACCTAGAGT.fastq.gz',
    'H3H7YDRXY_2#144456_ACTCGCTACGTCTAAT.fastq.gz'
]


class Test_align_fastq_with_plate_positions:

    def test_with_external_file(self):
        """
        Test align_fastq_with_plate_positions with external data.
        """
        index_sets = "methlab/data/nordborg_nextera_index_sets.csv"
        x = ml.align_fastq_with_plate_positions(file_list, index_sets, "test")
        assert isinstance(x, pd.DataFrame)
        assert x.shape[0] == 3
        assert all( x.keys() == ['sample', 'fastq_1', 'fastq_2'])
    
    def test_with_internal_file(self):
        """
        Test align_fastq_with_plate_positions with internal data.
        """
        x = ml.align_fastq_with_plate_positions(file_list, 'nordborg', "test")
        assert isinstance(x, pd.DataFrame)
        assert x.shape[0] == 3
        assert all( x.keys() == ['sample', 'fastq_1', 'fastq_2'])

    def test_correct_order(self):
        """
        confirm that the output gives files in the right order, even if the input does not
        """
        wonky_file_list = [file_list[i] for i in [3,4,5,0,1,2]]
        x = ml.align_fastq_with_plate_positions(
            wonky_file_list,
            'nordborg',
            "test"
            )
        fastq1 = x['fastq_1'].tolist()
        fastq2 = x['fastq_2'].tolist()

        assert all(
            [ fastq1[i][10] < fastq2[i][10] for i in [0,1,2] ]
        )
