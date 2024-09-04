import numpy as np
from warnings import warn
import numpy.ma as ma
import pandas as pd
import plotly.express as px
import allel
import h5py

def load_genotypes(input, reference, sample_name):
    """
    Import and merge test and reference data files.

    Import genotype data for one or more input samples and a panel of reference samples
    to compare to. Subset each so that the markers are really identical. Merge the
    arrays of genotype calls so that the data for the input appear first on the
    first axis of the genotype call arrays.

    Parameters
    ==========
    input: str
        Path to a VCF file containing genotype data for one or more samples to check
    reference: str
        Path to a HDF5 file containing genotype data for a panel of reference individuals
        to compare the input indivual against.
    sample_name: str
        Sample name for the individual to check. This must be present in the samples
        in the input VCF.

    Return
    ======
    An object of class geneticDistance.
    """
    # Read in the data files
    input_vcf = allel.read_vcf(input, samples=[sample_name])
    ref_hdf5  = h5py.File(reference, mode="r")

    ref_str_data = {
        'samples' : [ x.decode('utf-8') for x in ref_hdf5['samples'][:] ],
        'chr'     : [ x.decode('utf-8') for x in ref_hdf5['variants/CHROM'][:] ]
    }

    if sample_name not in input_vcf['samples']:
        raise ValueError("The sample name is not in the list of samples in the input VCF file.")
    else: 
        # Find the position of the individual to test
        sample_ix = np.where(input_vcf['samples'] == sample_name)[0][0]
        # Join vectors of sample names, with the test individual first
        new_samples = np.append(
            input_vcf['samples'][None,sample_ix],
            ref_str_data['samples']
            )        

    # Check that contig labels match
    chr_labels = {
        'input' : np.unique(input_vcf['variants/CHROM']),
        'ref'   : np.unique(ref_str_data['chr'])
    }
    if len(chr_labels['input']) != len(chr_labels['ref']):
        raise ValueError(
            "The number of unique contig labels do not match: the input VCF has {}, but the reference panel has {}.".
            format( chr_labels['input'], chr_labels['ref'] )
        )
    elif any( chr_labels['input'] != chr_labels['ref'] ):
        raise ValueError(
            "Contig labels do not match between the input and reference files."
        )
    
    # Make sure we only compare SNPs that are found in both datasets.
    # Concatenate chromosome labels and SNP positions
    snp_names = {
        'input' : [ str(chr) + ":" + str(pos) for chr,pos in zip(input_vcf['variants/CHROM'], input_vcf['variants/POS']) ],
        'ref'   : [ str(chr) + ":" + str(pos) for chr,pos in zip(ref_str_data['chr'],   ref_hdf5['variants/POS'][:]) ]
    }
    # Find the SNP position names that are common to both datasets
    matching_SNPs_in_both_files = np.intersect1d(
        snp_names['input'],
        snp_names['ref']
        )
    which_SNPs_to_keep = {
        "input" : [ x in matching_SNPs_in_both_files for x in snp_names['input'] ],
        "ref"   : [ x in matching_SNPs_in_both_files for x in snp_names['ref'] ]
    }


    # Append the genotype data for the test individual to the array of the reference panel
    new_geno = np.concatenate(
        (input_vcf['calldata/GT'][which_SNPs_to_keep['input'], sample_ix][:, np.newaxis],
        ref_hdf5['calldata/GT'][which_SNPs_to_keep['ref']]),
        axis=1
        )
    
    # Define an output before closing the Hdf5 file
    output = geneticDistance(
        samples = new_samples,
        chr = np.array(ref_str_data['chr'])[np.where(which_SNPs_to_keep['ref'])[0]],
        pos = ref_hdf5['variants/POS'][:][which_SNPs_to_keep['ref']],
        geno = new_geno
    )

    ref_hdf5.close()

    return output

class geneticDistance(object):
    """
    A simple class to compare genotype data genetic distances between individuals.  

    Parameters
    ==========
    samples: array
        Vector of length m giving names for each sample.
    chr: array
        Vector of length n giving chromosome labels for each SNP.
    pos: array
        Vector of length n giving SNP positions. Note that SNP positions are inherited from 
        skikit allel and give row numbers from the input VCF file rather than
        base-pair positions on each chromosome.
    geno: array
        m x n x 2 array of genotype data where axes index SNPs, individuals, and 
        homologous chromosomes.

    Attributes
    ==========
    samples: array
        Vector of M sample names. The first sample is the input individual to be
        compared to the remaining reference individuals.
    chr: array
        Vector of chromosome labels. These are imported from the reference panel.
    pos: array
        Vector of N SNP positions. These are imported from the reference panel.
    geno: array
        NxMx2 array of genotype data, where N is the number of SNPs and M is the
        number of samples.

    Methods
    =======
    split_into_windows
        Split a geneticDistance object into windows.
    pairwise_distance
        Calculate pairwise genetic distance between an input individual and all 
        reference individuals.
    
    """
    def __init__(self, samples, chr, pos, geno):
        self.samples = samples
        self.chr = chr
        self.pos = pos
        self.geno = geno

    def split_into_windows(self, window_size: int):
        """
        Split a geneticDistance object into windows.

        Splits the geneticDistance object into chromosomes, then into windows on each
        chromosome. It returns a dictionary of geneticDistance objects for each window.

        Parameters
        ==========
        window_size: int
            Window size in base pairs.

        Returns
        =======
        A dictionary of geneticDistance objects with an element for each window.
        Indexes are in the form "Chr:start-stop".
        """
        # Empty dict to store the output
        list_of_distance_objects = {}

        for chr in np.unique(self.chr):
            # Boolean array indexing items in this chromosome
            chr_ix = self.chr == chr
            
            # Array of starting positions for each window. 
            start_positions = np.arange(0, self.pos[chr_ix].max(), window_size)
            for start in start_positions:
                stop  = start + window_size
                # Index positions of SNPs within the current window
                window_ix = (self.pos[chr_ix] >= start) & (self.pos[chr_ix] < stop)
                # Create an object for each window.
                window_name = str(chr) + ":" + str(start) + "-" + str(stop)
                list_of_distance_objects[window_name] = geneticDistance(
                        samples = self.samples,
                        chr  = self.chr[chr_ix][window_ix],
                        pos  = self.pos[chr_ix][window_ix],
                        geno = self.geno[chr_ix][window_ix]
                    )
        
        return list_of_distance_objects
    
    def pairwise_distance(self, warn_about_missing_data=False):
        """
        Calculate pairwise genetic distance between an input individual and all 
        reference individuals.

        The input individual is always the first in the list of samples. Genetic
        distance is the number of allelic differences at each locus between each
        pair, summed over all loci. The calculation is done using masked arrays to
        account for missing data.

        Returns
        =======
        Vector of distances

        """
        masked_geno = ma.masked_array(self.geno, self.geno < 0)

        # Calculate differences at each locus
        per_locus_difference = abs(masked_geno.sum(2)[:,[0]] - masked_geno.sum(2)[:,1:]) / 2
        # Average over loci
        dxy = per_locus_difference.mean(0)

        if warn_about_missing_data and any(dxy.mask):
            warn("""

    Pairwise distance could not be calculated for one or more comparisons,
    probably because there is missing data at all SNPs.
    The following samples in the reference panel are affected:
    {}

    """.format(self.samples[1:][dxy.mask])
            )
            # Return a vector of -9 t indicate missing data
            return np.zeros(len(self.samples)-1) -9
        
        return dxy.data

def ibd_table(input:str, reference:str, sample_name:str, window_size:int):
    """
    Compare allele sharing across the genome.

    Calculate genetic distance between a test individual and a panel of
    reference genomes.

    Parameters
    ==========
    input: str
        Path to a VCF file containing genotype data for one or more samples to
        test
    reference: str
        Path to an HDF5 file containing genotype data for a panel of reference
        individuals to compare the test individual against.
    sample_name: str
        Sample name for the individual to check.
        This must be present in the samples in the input VCF.
    window_size: int
        Window size in base pairs.

    Returns
    =======
    DataFrame with a row for each window in the genome and a column for each 
    sample in the reference panel. Elements show genetic distance between the 
    test individual and each reference individual in a single window.
    """
    genetic_distance = load_genotypes(
            input = input,
            reference = reference,
            sample_name = sample_name
        )
    # Divide the genome into windows
    distances_in_windows = genetic_distance.split_into_windows(window_size)

    # Dataframe with a row for each window across the genome and a column for each sample in the reference panel.
    distance_array = pd.DataFrame(
        [ v.pairwise_distance() for v in distances_in_windows.values() ],
        columns = genetic_distance.samples[1:]
    )
    distance_array.insert(0, 'window', distances_in_windows.keys())

    return distance_array

def ibd_scores(ibd_table):
    """
    Mean IBD across the genome.

    Calculate mean genetic distance from a test individual to each of a panel of
    reference samples, ignoring windows where there was only missing data.

    Parameters
    ==========
    ibd_table: pd.DataFrame
        DataFrame with a row for each window in the genome and a column for each 
        sample in the reference panel. Elements show genetic distance between the 
        test individual and each reference individual in a single window.
        This is generated by ibd_table().

    Returns
    =======
    A DataFrame with a row for each candidate in the reference panel, and a 
    column indicating mean genetic distance over windows across the genome.
    Values closer to zero indicate that the sample is more likely to be a match.
    """
    # Coerce missing data to NaN for correct column means.
    ibd_table = ibd_table.replace(-9,np.NaN)

    # Get column-mean IBD for each candidate, allowing for missing data
    ibd_scores_for_each_candidate = np.array(
        [ np.nanmean(ibd_table[col]) for col in ibd_table.keys()[1:] ]
    )
    scores = pd.DataFrame({
        'candidate': ibd_table.keys()[1:],
        'score' : ibd_scores_for_each_candidate
    })
    return scores

def plot_ibd_table(ibd_table:pd.DataFrame, sample_name:str, expected_match:list=[], max_to_plot=20):
    """
    Plot allele sharing across the genome.


    Create a interactive line graph showing genetic distance from a test
    individual to each sample in a panel of reference individuals.

    Parameters
    ==========
    ibd_table: pd.DataFrame
        DataFrame with a row for each window in the genome and a column for each 
        sample in the reference panel. Elements show genetic distance between the 
        test individual and each reference individual in a single window.
        This is generated by ibd_table().
    sample_name: str
        Sample name for the individual to check.
        This must be present in the samples in the input VCF.
    expected_match: list
        List of sample names in the reference panel that are expected to be
        ancestors of the test individual.

    Returns
    =======
    Plotly figure object with subplots for each chromosome, showing window
    position along the x-axis and genetic distance from the test individual to
    each reference sample on the y-axis. Line colour indicates whether a sample
    is an expected parent or not. Rolling over the lines shows which sample is
    which.
    """

    
    # Coerce missing data to NaN for correct column means.
    ibd_table = ibd_table.replace(-9,np.NaN)

    # Identify the candidate names *not* among the top `max_to_plot` columns and remove
    # If `max_to_plot` is less than the number of candidates.
    if max_to_plot < ibd_table.shape[1]-1:
        # Get column-mean IBD for each candidate, allowing for missing data
        ibd_scores_for_each_candidate = np.array(
            [ np.nanmean(ibd_table[col]) for col in ibd_table.keys()[1:] ]
        )
        # Identify the candidate names *not* among the top `max_to_plot` columns
        ix = np.argpartition(ibd_scores_for_each_candidate, max_to_plot)[max_to_plot:] # index positions
        columns_to_drop = ibd_table.keys()[ix+1].to_list() # candidate names
        ibd_table = ibd_table.drop(columns=columns_to_drop) # drop the candidates

    # Make the table long
    ibd_table = ibd_table.melt(id_vars=['window'], var_name='candidate', value_name='distance')
    # Column indicating which candidates should be plotted a different colour.
    ibd_table['expected'] = ibd_table['candidate'].isin(expected_match)
    # Split the 'window' column up into separate columns for chromosome, start and stop positions
    ibd_table[['chr', 'window']] = ibd_table['window'].str.split(":", expand=True)
    ibd_table[['start', 'stop']] = ibd_table['window'].str.split("-", expand=True)
    # start and stop positions should be integers for sensible plotting.
    ibd_table['start'] = ibd_table['start'].astype(int)
    ibd_table['stop'] = ibd_table['stop'].astype(int)
    ibd_table['midpoint'] = (ibd_table['start'] + ibd_table['stop']) / 2

    fig = px.line(
        ibd_table,
        x="midpoint", y="distance", color="expected",
        title=sample_name,
        labels={
            'midpoint' : 'Position (bp)',
            'distance' : 'Genetic distance'        
        },
        hover_data=['candidate'],
        category_orders={'expected': [False, True]},
        facet_row = "chr"
        )
    fig.update_traces(mode="markers+lines")

    return fig