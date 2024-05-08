import numpy as np
from warnings import warn
import numpy.ma as ma
import allel
allel.__version__

"""
In:
    VCF file for input individuals to test
    Sample ID for the single individual to test
    VCF file for reference individuals
    Integer number of SNPs to look at. Have to round this
Out:
    CSV file with a row for each window and a column for each reference individual
        giving the genetic distance from the input sample to each reference.

Function to import two VCF files
    Subset input VCF to input sample only.
    Check that the two contain the same number of SNPs
    Concatenate genotype arrays?
Calculate distance from an input individual to each individual in the reference
"""

def import_vcfdistance(input, reference, sample_name):
    """
    Import and merge data from VCF files.

    Import VCF files for one or more input samples and a panel of reference samples
    to compare to. Subset each so that the markers are really identical. Merge the
    arrays of genotype calls so that the data for the input appear first on the
    first axis of the genotype call arrays.

    Parameters
    ==========
    input: str
        Path to a VCF file containing genotype data for one or more samples to check
    reference: str
        Path to a VCF file containing genotype data for a panel of reference individuals
        to compare the input indivual against.
    sample_name: str
        Sample name for the individual to check. This must be present in the samples
        in the input VCF.

    Return
    ======
    An object of class VcfDistance.
    """
    # Read in the VCF files
    input_vcf = allel.read_vcf(input)
    ref_vcf   = allel.read_vcf(reference)

    if sample_name in input_vcf['samples']:
        # Find the position of the individual to test
        sample_ix = np.where(input_vcf['samples'] == sample_name)[0][0]
        # Join vectors of sample names, with the test individual first
        new_samples = np.append(
            input_vcf['samples'][None,sample_ix],
            ref_vcf['samples']
            )
    else:
        raise ValueError("The sample name is not in the list of samples in the input VCF file.")

    # Check that contig labels match
    chr_labels = {
        'input' : np.unique(input_vcf['variants/CHROM']),
        'ref'   : np.unique(ref_vcf['variants/CHROM'])
    }
    if len(chr_labels['input']) != len(chr_labels['ref']):
        raise ValueError(
            "The number of unique contig labels do not match: the input VCF has {}, but the reference VCF has {}.".
            format( chr_labels['input'], chr_labels['ref'] )
        )
    elif any( chr_labels['input'] != chr_labels['ref'] ):
        raise ValueError(
            "Contig labels do not match between the input and reference VCF files."
        )
    
    # Make sure we only compare SNPs that are found in both VCF files.
    # Concatenate chromosome labels and SNP positions
    snp_names = {
        'input' : [ str(chr) + ":" + str(pos) for chr,pos in zip(input_vcf['variants/CHROM'], input_vcf['variants/POS']) ],
        'ref'   : [ str(chr) + ":" + str(pos) for chr,pos in zip(ref_vcf['variants/CHROM'],   ref_vcf['variants/POS']) ]
    }
    # Find the SNP position names that are common to both VCF files
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
        ref_vcf['calldata/GT'][which_SNPs_to_keep['ref']]),
        axis=1
        )

    return VcfDistance(
        samples = new_samples,
        chr = ref_vcf['variants/CHROM'][which_SNPs_to_keep['ref']],
        pos = ref_vcf['variants/POS'][which_SNPs_to_keep['ref']],
        geno = new_geno
    )






class VcfDistance(object):
    """
    A simple class to compare genotype data genetic distances between individuals
    from a VCF file.  

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
        Split a VcfDistance object into windows.
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
        Split a VcfDistance object into windows.

        Splits the VcfDistance object into chromosomes, then into windows on each
        chromosome. It returns a dictionary of VcfDistance objects for each window.

        Parameters
        ==========
        window_size: int
            Window size in base pairs.

        Returns
        =======
        A dictionary of VcfDistance objects with an element for each window.
        Indexes are in the form "Chr:start-stop".
        """
        # Empty dict to store the output
        list_of_vcfdistance_objects = {}

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
                list_of_vcfdistance_objects[window_name] = VcfDistance(
                        samples = self.samples,
                        chr  = self.chr[chr_ix][window_ix],
                        pos  = self.pos[chr_ix][window_ix],
                        geno = self.geno[chr_ix][window_ix]
                    )
        
        return list_of_vcfdistance_objects
    
    def pairwise_distance(self):
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

        if any(dxy.mask):
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
