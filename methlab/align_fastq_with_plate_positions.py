import pandas as pd
import re
import os
from warnings import warn

def align_fastq_with_plate_positions(
        mate1:list,
        mate2:list,
        sample_sheet:pd.DataFrame
        ):
    """
    Align raw sequence files with plate positions

    align_fastq_with_plate_positions looks up the adapter sequence in a lists of
    (paired) raw sequence files (usually .fastq.gz) in a data frame of adapter
    sequences to determine the position (row/column) of the sample in the 
    sequencing plate.

    This relies on grepping a string of only A, T, C or G of at least 8 letters
    from the file name of fastq files (basename only; the longer path is ignored).
    If such a string exists in the file name but is not an index, expect problems.

    Parameters
    ==========
    mate1: list
        List of paths to raw sequence files (usually .fastq.gz) for the first
        mate pairs. These should all be from a single sequencing plate, or else
        there will be multiple matches to each row/column position.
    mate2: list
        As for mate1, but for second mate pairs. Matching pairs do not need to
        be in the same order, but there should be one and only one file name in
        mate2 with the same index as mate2.
    sample_sheet: DataFrame
        Pandas dataframe with a row for each sample giving information about
        each sample. At a minimum this must contain columns 'index1' and 'index2'
        giving the forwards and reverse indices to look up.

    Returns
    =======
    The original data frame with additional columns 'fastq1' and 'fastq2' giving
    paths to the fastq files.
    """
    # Check mate1 and mate2 are the same length
    if len(mate1) != len(mate2):
        raise ValueError("mate1 and mate2 are not the same length.")
    # Check column headers
    check_col_names = all([col_name in sample_sheet.keys() for col_name in ['index1', 'index2'] ])
    if not check_col_names:
        raise ValueError("`sample_sheet` should contain at least the headers 'index1' and 'index2'")

    # Merge adapters into a single sequence.
    sample_sheet['seq_combined'] = sample_sheet['index1'] + sample_sheet['index2']

    # Dataframe with a row for each input file giving the two indices pasted
    # together, with both full paths
    input_files_df = pd.DataFrame({
        # Pull the indices from the 
        'seq_combined' : [re.findall('[ACTG]{8,}', os.path.basename(path_name))[0] for path_name in mate1],
        'fastq1' : mate1,
        'fastq2' : mate2
    })

    duplicate_indices = input_files_df['seq_combined'].duplicated()
    if any(duplicate_indices):
        warn(
            UserWarning("One or more pairs of input files have duplicated indices.")
        )

    sample_sheet = sample_sheet.merge(input_files_df, how='left', on='seq_combined')

    return sample_sheet