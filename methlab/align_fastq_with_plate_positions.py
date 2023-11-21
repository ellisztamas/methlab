import pandas as pd
import re
import numpy as np

import sys

if sys.version_info < (3, 9):
    # importlib.resources either doesn't exist or lacks the files()
    # function, so use the PyPI version:
    import importlib_resources
else:
    # importlib.resources has files(), so use that:
    import importlib.resources as importlib_resources

def import_index_set(index_set:str):
    """
    Import a set of adapter indices

    Import a data frame giving positions and sequences of a set of adapater
    indices. These can be from a set of files stored internally in the 
    methlab package, or from a path to an external file.

    The full list of adapters at the VBC NGS facility is given on the forskalle
    portal: https://ngs.vbcf.ac.at/forskalle3/account/adaptors
    These sets are implemented internally:
    - 'nordborg': Nordborg Nextera INDEX set 1-8
    - 'dualxt' : Unique Nextera Dual XT. Note: this needs to be confirmed by
        checking it returns the correct answer for a known plate - it is
        possible that one or more index sets are the wrong way round.

    For other sets the exact relationship between adapter sequence and position
    is left as a cunning riddle to the user.

    Parameters
    ==========
    index_set: str
        Name of the index set. Provide a name of an internal index set, or a
        path to an external CSV file. Internal sets may be 'dualxt' or
        'nordborg'; see above. External files should contain columns called
        'row', 'col', 'seq1' and 'seq2'.

    Returns
    =======
    pd.DataFrame
        A data frame giving row and column position, sequences of each pair of
        adapters, and potentially additional columns.
    """
    package_index_sets = {
        'nordborg' : "nordborg_nextera_index_sets.csv",
        'dualxt'   : "unique_nextera_dual_xt.csv"
    }

    if index_set in ['nordborg', 'dualxt']:
        pkg = importlib_resources.files("methlab")
        adapter_indices = pd.read_csv(pkg/"data"/package_index_sets[index_set])
    else:
        adapter_indices = pd.read_csv(index_set)
    
    # Check column headers
    check_col_names = all([col_name in adapter_indices.keys() for col_name in ['row', 'col', 'seq1', 'seq2'] ])
    if not check_col_names:
        raise ValueError("`adapter_indices` should contain at least the headers 'row', 'col', 'seq1' and 'seq2'")
    
    # Merge adapters into a single sequence.
    adapter_indices['seq_combined'] = adapter_indices['seq1'] + adapter_indices['seq2']

    return adapter_indices

def align_fastq_with_plate_positions(input_files:list, adapter_indices:str, prefix:str=""):
    """
    Align raw sequence files with plate positions

    align_fastq_with_plate_positions looks up the adapter sequence in a list of
    raw sequence files (usually .fastq.gz) in a data frame of adapter sequences
    to determine the position (row/column) of the sample in the sequencing plate.

    Parameters
    ==========
    input_files: list
        List of paths to raw sequence files (usually .fastq.gz). These should
        all be from a single sequencing plate, or else there will be multiple
        matches to each row/column position.
    adapter_indices: str
        Name of the index set. Provide a name of an internal index set, or a
        path to an external CSV file. Two index sets are
        available internally in the package and can be used by passing
        'nordborg' for the Nordborg Nextera INDEX set 1-8 or 'dualxt' for the
        Unique Nextera Dual XT indices. External files should contain columns
        called 'row', 'col', 'seq1' and 'seq2'. If there are multiple index sets
        (for example there are 8 sets for the Nextera sets), give a row for all 
        possible sets. 
    prefix: str
        Prefix name to be given to the samples. This is usually a label for the 
        plate name. This is appended to the position label in the output.

    Returns
    =======
    pd.DataFrame
        A dataframe giving sample name, path to mate pair 1 and path to mate pair 2.
        If data a single-end, the path to mate pair 2 will be blank.
    """
    adapter_indices = import_index_set(adapter_indices)

    # For each fastq file, find the position of the matching adapter sequence in the adapter indices.
    ix = []
    for path_name in input_files:
        input_adapter_sequence = re.findall('[ACTG]+', path_name)[0]
        row_number = np.where(
            adapter_indices['seq_combined'].str.match(input_adapter_sequence)
            )[0][0]
        ix = ix + [row_number]
    # Reorder adpater_indices to match the order of filenames, and add the file_names
    updated_adapter_indices = adapter_indices.loc[ix]
    updated_adapter_indices['file_name'] = input_files
    
    # Single column giving row/column position
    updated_adapter_indices['pos'] = updated_adapter_indices['row'].astype(str) + updated_adapter_indices['col'].astype(str)
    # list of unique positions
    position_labels = updated_adapter_indices['pos'].unique()
    
    # Create a dataframe giving desired sample name, plus names of paired fastq files
    # If data are unpaired, the values for file 2 are left blank.
    sample_sheet = []
    for label in position_labels:
        # For each plate position, return a list of all matching filenames
        matching_filenames = updated_adapter_indices.loc[updated_adapter_indices['pos'] == label]['file_name']
        matching_filenames = sorted( matching_filenames.to_list() )
        # There will usually be two files corresponding to paired reads.
        # If there is only one matching filename, add a blank element for a fictional read pair
        if len(matching_filenames) == 1:
            matching_filenames = matching_filenames + ['']
        # If there are more than two, raise and exception
        elif len(matching_filenames) > 2:
            raise ValueError("There are more than two samples matching position " + label)
        
        # Create a row of a data frame giving the name of the sample, plus all filenames
        sample_name = [prefix + label]
        sample_sheet = sample_sheet + [sample_name + matching_filenames]

    # Export as a data frame.
    sample_sheet = pd.DataFrame(sample_sheet, columns = ['sample','fastq_1','fastq_2'])
    
    return sample_sheet.sort_values('sample')