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
    if len(mate1) == 0:
        raise ValueError("List of fastq files is empty!")
    # Check column headers
    check_col_names = all([col_name in sample_sheet.keys() for col_name in ['index1', 'index2'] ])
    if not check_col_names:
        raise ValueError("`sample_sheet` should contain at least the headers 'index1' and 'index2'")

    # Merge adapters into a single sequence.
    sample_sheet['combined_indices'] = sample_sheet.loc[:,'index1'] + sample_sheet.loc[:,'index2']

    # Dataframes giving indices and paths to input files for mates 1 and 2.
    # It is assumed indices are at least 8 nucleotides
    indices_mate1 = pd.DataFrame({
        # Pull the indices from the 
        'combined_indices' : [re.findall('[ACTG]{8,}', os.path.basename(path_name))[0] for path_name in mate1],
        'fastq1' : mate1
    })
    indices_mate2 = pd.DataFrame({
        # Pull the indices from the 
        'combined_indices' : [re.findall('[ACTG]{8,}', os.path.basename(path_name))[0] for path_name in mate2],
        'fastq2' : mate2
    })
    # Merge into a single dataframe with an inner join
    merged_indices_fastq = pd.merge(indices_mate1, indices_mate2, how='inner', on='combined_indices')
    # Because this is an inner join, if one or more indices do not match between 
    # mate1 and mate2 then merged_indices_fastq will have fewer rows than the number
    # of files. Check this.
    if merged_indices_fastq.shape[0] != len(mate1):
        raise ValueError("One or more indices is present in mate1 but not mate2, or vice versa")
    
    sample_sheet = sample_sheet.merge(merged_indices_fastq, how='left', on='combined_indices')

    return sample_sheet