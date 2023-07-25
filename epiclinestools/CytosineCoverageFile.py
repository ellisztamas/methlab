import pandas as pd


class CytosineCoverageFile(object):
    """
    A class to work with Bismark coverage2cytosine files

    Parameters
    ==========
    path: str
        Path to a cytosine2coverage-format file from Bismark. This should have
        no header (this will be added), but seven columns giving chromosome,
        base-pair of each cytosine position on the chromosome, strand (+/-), 
        number of methylated reads, number of unmethylated reads, sequence
        context (CG, CHG, CHH) and tricnucleotide context (for example, for CHG
        methylation this could be CAG, CCG or CTG).
    
    Examples
    ========
    # Example coverage file
    path="tests/test_data/test_coverage2cytosine.txt.gz"

    # Load the file into memory
    c2c = CytosineCoverageFile(path)

    # Subset 10 bp on Chr1
    chr = "Chr1"
    start = 3000
    stop = 3010
    sub = c2c.subset(chr, start, stop)




    """
    def __init__(self, path:str):
        self.path = path
        self.file = self.read()

    def read(self):
        """
        Read in a cytosine2coverage file.

        Load a cytosine2coverage file into memory from the path given in the 
        class instance. If it is compressed (which it probably is) there is no
        need to unzip it; pandas handles this automatically.
        """
        file = pd.read_csv(
            self.path, sep="\t",
            names= ['chr', 'pos', 'strand', 'meth', 'unmethylated', 'context', 'trinucleotide'],
            dtype={
                'chr' : 'category',
                'pos' : 'Int64',
                'strand': "category",
                'methylated' : 'Int64',
                'unmethylated': 'Int64',
                'context': 'category',
                'trinucleotide' : "category"
            }
            )
        return file

    def subset(self, chr, start:int, stop:int):
        """
        Subset the file based on chromosome, start and stop position.

        Parameters
        ==========
        chr: str or int
            Label for the chromosome label.
        start: int
            Starting base position
        stop: stop
            Stop base position
        
        Returns
        =======
        Pandas dataframe which is a row subset of self.file.
        """
        return self.file.loc[
            (self.file['chr'] == chr) & 
            (self.file['pos'] >= start) & 
            (self.file['pos'] <= stop)
            ]

    def count_reads(self, data):
        """
        Helper function to count methylated and unmethylated reads.

        Parameters
        ==========
        data: pd.DataFrame
            A whole or partial coverage file containing at least column headers
            'context', 'meth' and 'unmethylated'.

        Returns
        =======
        Pandas dataframe specifying sequence context, number of methylated and
        unmethylated reads, and number of cytosines.
        """
        # Read counts and cytosine number in each context
        mC_read_counts = data.groupby('context').sum(numeric_only = True)[['meth', 'unmethylated']]
        mC_read_counts['ncytosines'] = data.groupby('context').size()
        # Add a row giving totals
        mC_read_counts.loc['total',:] = mC_read_counts.sum(0).tolist()
        mC_read_counts = mC_read_counts.astype(int)

        return mC_read_counts
    
    def methylation_over_features(self, chr, start, stop, names= None):
        """
        Methylated reads across annotated features

        Counts the number of methylated and unmethylated reads and number of 
        cytosines in CG, CHG and CHH contexts in each of multiple annotated
        features (e.g. genes, TEs).

        Note that values in chr must all match a chromosome in the coverage file.

        Parameters
        ==========
        chr: vector
            Vector of chromosome labels to look up in the coverage file.
        start: vector
            Vector of positions indexing the base-pair position of the start of
            each feature.
        stop: vector
            Vector of positions indexing the base-pair position of the end of
            each feature.
        names: vector, optional
            Vector of names for each feature.

        Returns
        =======
        Pandas dataframe specifying feature, sequence context, number of 
        methylated and unmethylated reads, and number of cytosines.

        Example
        =======
        # Example annotation file using the first ten lines of the TAIR10 annotation
        gff_file = pd.read_csv(
            "tests/test_data/test_TAIR10_GFF3_genes_transposons.gff",
            sep="\t",
            names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
            ).iloc[1:9] # Skip the first row, because it defines the whole chromosome
        # Example coverage file
        path="tests/test_data/test_coverage2cytosine.txt.gz"
        c2c = CytosineCoverageFile(path)
        
        meth_counts = c2c.methylation_over_features(
            chr = gff_file['seqid'],
            start = gff_file['start'],
            stop = gff_file['end']
            )    
        """
        check_chr_labels = any(~chr.isin(self.file['chr']))
        if check_chr_labels:     
            raise ValueError("One or more values in `chr` is not found in the chromosome labels in the coverage file")

        if names is None:
            names = ['feature' + str(x) for x in range(len(chr))]
        
        coords = zip(names, chr, start, stop)

        features = {}
        for name, chr, start, stop in coords:
            this_feature = self.subset(chr, start, stop)
            features[name] = self.count_reads(this_feature)

        features = pd.concat(features).reset_index()
        features = features.rename(columns={'level_0': 'id'})

        return features

    def conversion_rate(self, return_proportion = True):
        """
        Calculate methylation on each chromosome. Conversion rate is 1 minus 
        these proportions.

        Parameters
        ==========
        chr: str, probably
            Probably a string indicating which chromosome to calculate 
            conversion rate on. Defaults to the chloroplast.
        
        Returns
        =======
        Pandas dataframe specifying chromosome, sequence context, proportion of
        methylated and unmethylated cytosines, and number of cytosines. If 
        `return_counts` is `True`, read counts are returned instead of
        proportions.

        Example
        =======
        # Example coverage file
        path="tests/test_data/test_coverage2cytosine.txt.gz"
        c2c = CytosineCoverageFile(path)

        # Return proportion of methylation on each chromosome
        c2c.conversion_rate()

        # Disable the conversion to proportions and return raw read counts
        c2c.conversion_rate(return_proportion = False)
        """
        chr_labels = self.file['chr'].unique()

        chromosomes = {}
        for chr in chr_labels:
            this_chr = self.file.loc[self.file['chr'] == chr]
            chromosomes[chr] = self.count_reads(this_chr)

        chromosomes = pd.concat(chromosomes).reset_index()
        chromosomes = chromosomes.rename(columns={'level_0': 'id'})

        # Convert to proportion methylated
        if return_proportion:
            chromosomes['n_reads']      = chromosomes['meth'] + chromosomes['unmethylated']
            chromosomes['meth']         = chromosomes['meth'].astype(float)         / chromosomes['n_reads']
            chromosomes['unmethylated'] = chromosomes['unmethylated'].astype(float) / chromosomes['n_reads']

        return chromosomes
    
    def methylation_in_windows(self, window_size:int):
        """
        Count methylated reads in fixed windows

        Counts the number of methylated and unmethylated reads and number of 
        cytosines in CG, CHG and CHH contexts in windows of fixed size across the
        genome.

        To do: consider allowing for only a subset of chromosomes.

        Parameters
        ==========
        window_size: int
            Window size in base pairs

        Returns
        =======
        Pandas dataframe giving chromosome label, start position of the window,
        sequence context, number of methylated and unmethylated reads, and total
        number of cytosines.

        Example
        =======
        path="tests/test_data/test_coverage2cytosine.txt.gz"
        c2c = CytosineCoverageFile(path)

        # Return proportion of methylation on each chromosome
        c2c.methylation_in_windows(window_size)
        """
        chr_labels = self.file['chr'].unique()


        reads_over_chromosomes = {}

        for chr in chr_labels:
            this_chr = self.file.loc[self.file['chr'] == chr]
            reads_for_this_chr = {}
            for i in range(0, this_chr['pos'].max(), window_size):
                this_window = this_chr.loc[(this_chr['pos'] >= i) & (this_chr['pos'] < i+window_size)]
                reads_for_this_chr[i] = self.count_reads(this_window)
            reads_over_chromosomes[chr] = pd.concat(reads_for_this_chr)

        reads_over_chromosomes = pd.concat(reads_over_chromosomes).reset_index()
        reads_over_chromosomes = reads_over_chromosomes.rename(columns={'level_0': 'chr', 'level_1':'start'})

        return reads_over_chromosomes