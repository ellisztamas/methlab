==========================================
Tools to work with cytosine coverage files
==========================================

Bismark's ``coverage2cytosine`` tools creates very useful cytosine coverage
reports that give the coverage of methylated and unmethylated reads at each
cytosine. The class ``CytosineCoverageFile`` provides various tools for 
extractin information of interest from these files.

Import a coverage file
======================

To import a coverage file to memory give a path to a cytosine2coverage-format file from Bismark.
standard format is to have no header (this will be added), but seven columns giving
chromosome, base-pair of each cytosine position on the chromosome, strand (+/-), 
number of methylated reads, number of unmethylated reads, sequence context
(CG, CHG or CHH) and tricnucleotide context (for example, for CHG methylation
this could be CAG, CCG or CTG).

.. code-block:: python

    import methlab as ml    path="tests/test_data/test_coverage2cytosine.txt.gz"
    c2c = CytosineCoverageFile(path)
    # Access the whole file, if you really must (it's big!)
    c2c.file

Conversion rate on each chromosome
==================================

Calculate overall methylation on each chromosome separately.
This is especially useful for calculating conversion rates by quantifying 
methylation on the chloroplast

.. code-block:: python

    # Return proportion of methylation on each chromosome
    c2c.conversion_rate()

You can also return raw read counts instead of proportions:

.. code-block:: python

    c2c.conversion_rate(return_proportion = False)

Here is an example of the output for chromosome 1 only.
Notice that the format is the proportion of *methylated* rather than *unmethylated*
cytosines, so it isn't really conversion rate in that the strict sense.

.. parsed-literal::

    id context      meth  unmethylated  ncytosines  n_reads
    Chr1      CG  0.356164      0.643836          82       73
    Chr1     CHG  0.333333      0.666667          68       54
    Chr1     CHH  0.144715      0.855285         818      615
    Chr1   total  0.179245      0.820755         968      742
    
Methylation in windows
======================

Count the number of methylated and unmethylated reads and number of  cytosines
in CG, CHG and CHH contexts in windows of fixed size across the genome.

.. code-block:: python

    # Methylation in 150-bp windows
    windows = c2c.methylation_in_windows(150)

Here is an example of the output. You can see that windows are indxed by their
start position.

.. parsed-literal::

        chr  start context  meth  unmethylated  ncytosines
    0   Chr1      0      CG    22             4          30
    1   Chr1      0     CHG    15             8          24
    2   Chr1      0     CHH    48           228         273
    3   Chr1      0   total    85           240         327
    4   Chr1   1000      CG     0             9          14
    5   Chr1   1000     CHG     2            17          28
    6   Chr1   1000     CHH    14           115         275
    7   Chr1   1000   total    16           141         317
    8   Chr1   2000      CG     4            34          38
    9   Chr1   2000     CHG     1            11          16
    10  Chr1   2000     CHH    27           183         270
    11  Chr1   2000   total    32           228         324

Methylation over annotated features
===================================

Count the number of methylated and unmethylated reads and number of cytosines
in CG, CHG and CHH contexts in each of multiple annotated features
(e.g. genes, TEs). Usually such information would be described by a file in 
``.bed`` or ``.gff`` format, but to keep it general we just use vectors of 
start and stop coordinates that are easily extracted from such files.

Here is an example using the first ten lines of the TAIR10 annotation file

.. code-block:: python
    
    gff_file = pd.read_csv(
        "tests/test_data/test_TAIR10_GFF3_genes_transposons.gff",
        sep="\t",
        names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        ).iloc[1:9] # Skip the first row, because it defines the whole chromosome
    # Example coverage file
    
    meth_counts = c2c.methylation_over_features(
        chr = gff_file['seqid'],
        start = gff_file['start'],
        stop = gff_file['end']
        )

Here is the output for the first three features.

.. parsed-literal::
            id context  meth  unmethylated  ncytosines
    0   feature0      CG     2            14          20
    1   feature0     CHG     0             4           8
    2   feature0     CHH     6            78         106
    3   feature0   total     8            96         134
    4   feature1      CG     2            14          20
    5   feature1     CHG     0             4           8
    6   feature1     CHH     6            78         106
    7   feature1   total     8            96         134
    8   feature2      CG     2            13          18
    9   feature2     CHG     0             4           6
    10  feature2     CHH     3            53          64
    11  feature2   total     5            70          88

This shows the generic default output for the feature name (the ``id`` column)
because the there are no useful IDs to use in the GFF file.
You can optionally supply a vector of names using the argument ``names``.
