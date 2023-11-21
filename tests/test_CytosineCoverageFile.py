import pandas as pd
import pytest

import methlab as ml
def test_read_CytosineCoverageFile():
    """
    Check reading in a file
    """
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)
    assert c2c.path == path
    assert c2c.file.shape == (10260, 7)

def test_subset():
    """
    Test subsetting self.file for a single feature using chromosome, start and 
    stop positions.
    """  
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)

    chr = "Chr1"
    start = 1000
    stop = 1010

    sub = c2c.subset(chr, start, stop)
    assert sub.shape == (2,7)

def test_methylation_over_features():
    """
    Test the function to calculate methylated and unmethylated reads over 
    features in an annotation file
    """

    # Example annotation file using the first ten lines of the TAIR10 annotation
    gff_file = pd.read_csv(
        "tests/test_data/test_TAIR10_GFF3_genes_transposons.gff",
        sep="\t",
        names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        ).iloc[1:9] # Skip the first row, because it defines the whole chromosome
    
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)
    
    meth_counts = c2c.methylation_over_features(
        chr = gff_file['seqid'],
        start = gff_file['start']-1000,
        stop = gff_file['end']-1000
        )
    # There should be 3 rows for every line in the GFF file
    meth_counts.shape[0] == 4 * gff_file.shape[0]
    meth_counts['id'].unique() == all(['feature' + str(i) for i in range(8)])

    # This shortens chromosome labels from "Chr1" to "1", which doesn't match
    # what is in the coverage file, and should raise an exception
    # with pytest.raises(Exception) as e_info:
    #     c2c.methylation_over_features(
    #         chr = gff_file['seqid'].str.slice(3,4),
    #         start = gff_file['start']-1000,
    #         stop = gff_file['end']-1000
    #         )

def test_conversion_rate():
    """
    Check this works over chromosomes
    
    To do: test with another coverage file that has more than one chromosome!
    """
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)
    cr = c2c.conversion_rate()
    assert all( cr['context'].isin(['CG', "CHG", "CHH", "total"]) )
    # Check it returns entries for all chromosomes
    chr_names = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrC', 'ChrM', 'Lambda_NEB', 'pUC19']
    assert all( cr['id'].isin(chr_names))

    # Check this works for a subset of chromosomes.
    cr_chloroplast = c2c.conversion_rate(chr_labels=['ChrC', "Lambda_NEB"])
    assert all( cr_chloroplast['id'].isin(['ChrC', "Lambda_NEB"]))

    # If a chromosome label is not in the cytosine file, return NaN
    # I don't know if this should return an error
    with pytest.warns(UserWarning):
        chr6 = c2c.conversion_rate(chr_labels=['Chr6'])

    # Commented out because pUC19 has zero reads at all and returns NaN
    # Need to decide how to handle that.
    # assert all( 
    #     cr['unconverted'] + cr['converted'] == 1
    #     )
    
    # Check return_proportion gives counts
    cr2 = c2c.conversion_rate(return_proportion=False)
    assert all(
        cr2.dtypes == ['object', 'object','int64', 'int64', 'int64']
    )
    # Commented out because pUC19 has zero reads at all and returns NaN
    # Need to decide how to handle that.
    # assert all(
    #     cr2['unconverted'] + cr2['converted'] > 1
    # )

def test_count_reads():
    """
    Check that count_reads does the totals properly.
    """
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)

    counts = c2c.count_reads(c2c.file)

    assert counts.iloc[0:3]['unconverted'].sum() == counts.iloc[3]['unconverted']
    assert counts.iloc[0:3]['converted'].sum() == counts.iloc[3]['converted']
    assert counts.iloc[0:3]['ncytosines'].sum() == counts.iloc[3]['ncytosines']

def test_methylation_in_windows():
    """
    Check basic functionalty of methylat_in_windows.
    """
    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = ml.CytosineCoverageFile(path)
    
    mc_windows = c2c.methylation_in_windows(1000)
    assert len(mc_windows['chr'].unique()) == 9

    chrC = c2c.methylation_in_windows(1000, chr_labels=['ChrC'])
    assert all(chrC['chr'] == "ChrC")