"""Console script for methlab."""
from align_fastq_with_plate_positions import align_fastq_with_plate_positions
import argparse
import sys

def align(input_files, adapter_indices, prefix=""):
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+', help="Array of fastq files to be aligned")
    parser.add_argument(
        'adapters', type=str,
        help="CSV files giving row and column positions, with the corresponding adaper sequence."
        )
    parser.add_argument(
        'prefix', type=str, required=False,
        help='Prefix name to be given to the samples. This is usually a label for the plate name. This is appended to the position label in the output.'
        )
    args = parser.parse_args()

    align_fastq_with_plate_positions(args.input_files, args.adapter_indices, args.prefix)

# def main():
#     """Console script for methlab."""
#     parser = argparse.ArgumentParser()
#     parser.add_argument('_', nargs='*')
#     args = parser.parse_args()

#     print("Arguments: " + str(args._))
#     print("Replace this message by putting your code into "
#           "methlab.cli.main")
#     return 0


# if __name__ == "__main__":
#     sys.exit(main())  # pragma: no cover
