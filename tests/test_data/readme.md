# Sample files for running unit tests

## VCF files

Sample VCF files for comparing genotyped samples against a reference panel:
```
reference_panel.vcf.gz
panel_to_test.vcf.gz
reference_panel_chr1.vcf.gz
```
These are taken from the methylation crosses project in
`/groups/nordborg/projects/epiclines/002.pedigree.
The reference panel is the four parents of the crosses performed (1158, 6024, 
6184, 8249) taken from Fernandos VCF.
The test panel are 10 samples of parents, F1s, F2s, F3s, a negative control and 
one positive control(Columbia).

The first two files should be the first 100kb on chromosomes 1 and 2.
The third is the same, but chromosome 1 only, to check that the program fails
if chromsome labels don't match.

Code to create the files:
```
bcftools view \
    --regions Chr1:1-100000,Chr2:1-100000 \
    --samples 1158,6024,6184,8249 \
    --output tests/test_data/reference_panel.vcf.gz \
    /groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_validate_genotyping/output/vcf/parents_only.vcf.gz

bcftools view \
    --regions Chr1:1-100000,Chr2:1-100000 \
    --samples 1158_2,6024_2,6184_1,8249_1,blank_2021-015_A4,Col0_2021-015_F9,F1.23.01,F1.26.03,F2.16.042,F3.03.001 \
    --output tests/test_data/panel_to_test.vcf.gz \
    /groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_validate_genotyping/output/vcf/pedigree_genotype_calls_against_TAIR10.vcf.gz

bcftools view \
    --regions Chr1:1-100000 \
    --samples 1158,6024,6184,8249 \
    --output tests/test_data/reference_panel_chr1.vcf.gz \
    /groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_validate_genotyping/output/vcf/parents_only.vcf.gz

```

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