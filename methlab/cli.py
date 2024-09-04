"""Console script for methlab."""
import argparse
import pandas as pd
from methlab.ibdpainting import *

# https://betterprogramming.pub/build-your-python-script-into-a-command-line-tool-f0817e7cebda

import argparse

def ibdpainting(input, reference, sample_name, window_size, expected_match, outdir, keep_ibd_table, max_to_plot):
    # Data frame of IBD at all positions across the genome, and the plot of this
    itable = ibd_table(input, reference, sample_name, window_size)
    scores = ibd_scores(itable)
    fig = plot_ibd_table(itable, sample_name, expected_match, max_to_plot)
    
    if keep_ibd_table:
        itable.to_csv(outdir + "/" + sample_name + "_ibd_table.csv", index=False)
    
    scores.to_csv(outdir + "/" + sample_name + "_ibd_scores.csv", index=False)
    fig.write_html(outdir + "/" + sample_name + "_plot_ibd.html")
    
  
def main():
    parser = argparse.ArgumentParser(description='methlab')

    # IBD painting
    subparsers = parser.add_subparsers(
        dest='command',
        help='Tools to automate tasks related to the analysis of DNA methylation variation.',
        required=True
        )
    parser_ibdpainting = subparsers.add_parser('ibdpainting',
        help="Calculate the genetic distance between a sample of interest and a panel of reference genotypes in windows across the genome"
        )
    parser_ibdpainting.add_argument('--input',
        help='Path to a VCF file containing genotype data for one or more samples to check.'
        )
    parser_ibdpainting.add_argument('--sample_name',
        help ='Sample name for the individual to check. This must be present in the samples in the input VCF.'
    )
    parser_ibdpainting.add_argument('--reference',
        help="Path to an HDF5 file containing genotype data for a panel of reference individuals to compare the input indivual against. This should be the output of allel.vcf_to_hdf5()"
    )
    parser_ibdpainting.add_argument('--window_size',
        type=int, default=20000,
        help="Integer window size in base pairs."
    )
    parser_ibdpainting.add_argument('--expected_match',
        help="Optional list of sample names in the reference panel that are expected to be ancestors of the test individual.",
        nargs = "+", required=False
    )
    parser_ibdpainting.add_argument('--outdir',
        help="Directory to save the output."
    )
    parser_ibdpainting.add_argument('--keep_ibd_table', 
        help="Optional ",
        action=argparse.BooleanOptionalAction
    )
    parser_ibdpainting.add_argument('--max_to_plot', 
        help="Optional number of the best matching candidates to plot so that the HTML files do not get too large and complicated. Ignored if this is more than the number of samples. Defaults to 20.",
        type=int, default = 20
    )
    parser_ibdpainting.set_defaults(func=ibdpainting)

    args = parser.parse_args()
    args_ = vars(args).copy()
    args_.pop('command', None)
    args_.pop('func', None)
    args.func(**args_)

if __name__ == '__main__':
    main()