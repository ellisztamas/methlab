# Sample files for running unit tests

## cytosine2coverage file

Sample file for testing the package is a subset of a real file with only ~10,000 lines:
```
test_coverage2cytosine.txt
```

Created by subsetting the following file as follows
```
# In python

path = "/groups/nordborg/projects/epiclines/006.quality_control/01_2022_bisulphite_protocol/04_output/30x_col0/bismark_meths/cx_report/220842_ATGTTGTTGGCAATCTATGA_S9_R1_001_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz"
c2c = pd.read_csv(
    path, sep="\t",
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

c2c.loc[c2c['pos'] < 3000].to_csv(
    "tests/test_data/test_coverage2cytosine.txt", index=False, sep="\t", header=False
    )

# Command line
gzip tests/test_data/test_coverage2cytosine.txt
```

## GFF annotation

Sample GFF file.

This is the first 10 lines of the TAIR10 annotation, created thus:

```
head /groups/nordborg/common_data/TAIR10/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff > /groups/nordborg/projects/epiclines/001.library/methlab/tests/test_data/test_TAIR10_GFF3_genes_transposons.gff
```

## gene_read_counts.csv

Example file for testing calling methylation state.
This is the output of CytosineCoverageFile.methylation_over_features applied to
the first ten genes from the TAIR10 annotation.