"""""""""""""""""""""""""""""""""
Align file names to 96-well plate
"""""""""""""""""""""""""""""""""

This page illustrates how to create a table listing sample names with paths to
their corresponding raw sequence data files for sequence data from the VBC
NGS facility.

This is a difficult problem to solve reliably because the formats keep changing,
so user discretion is essential.

Sequence file names
===================

When we submit 96-well plates for sequencing to the VBC NGS facility, we typically
submit something like an Excel sheet giving rows/columns of the plate and the
biological sample that should be in each well. What we get back from the NGS
facility is a mass of bam or fastq files that look something like this:

.. parsed-literal::

    343561_TTATACGCGAGCCGTCTGTT_S22_R1_001.fastq.gz

Older samples might look something like this:

.. parsed-literal::
  
    H3H7YDRXY_1\\#144456_ACTCGCTACGTCTAAT.bam

How is one meant to determine which sample in the original plate each file is
meant to correspond to?

What do the file names mean?
============================

Unfortunately, the naming system of files has changed through time, so they might
not always look like the ones above. In this case at least you can see:

* ``343561`` is the facility sample number. This usually indicates and ID for
  the plate. You can use this to track down the
  facility's data on the sequencing run (in this case, for example: 
  https://ngs.vbcf.ac.at/forskalle3/samples/343561).
  Confusingly, the facility also has a 'request number', which looks very similar.
* ``TTATACGCGAGCCGTCTGTT`` gives the **adapter index** sequence for this sample.
  This comprises two 8-or-more nucleotide sequences that together give a unique
  identifier for the row/column position in a 96-well plate. There may be
  multiple combinations for separate plates so that these can be run on a single
  flow cell.
* ``S22`` gives a number from 1 to 96 indicating the position in the well.
* ``R1`` indicates that these are paired reads, and that this is the file for
  mate-pair 1. Data for mate pair 2 would be found in a file with a similar name,
  but with ``R2`` in this position.
* ``001``: What does this mean? Only the Gods know.

On the older file:

* ``H3H7YDRXY`` is the flow cell on which sequences were run.
* ``_1`` is some kind of subset of the data on that flow cell.
  For example, 1 and 2 here might indicate either end of paired-end data.
  Note that sometimes your data might be combined with someone else's data, so
  you might have 3 and 4, or some other complicated combination of data.
  It's best to ask Almudena or Viktoria if you aren't sure.
* ``144456`` is the facility sample number.
* ``ACTCGCTACGTCTAAT`` gives the **adapter index** sequence for this sample.

As you can see, only the information the two formats have in common are the 
**sample ID** (plate) and adapter sequences, so we will use this information to
piece together which sample is which. Note that other formats may, or may yet
exist!

Which adapters do I have? 
=========================

For data created from 2023 or later you probably have the Nextera Dual XT
adapters, which are four sets of 96 unique pairs of adapters.
You can view them 
`here <https://docs.google.com/spreadsheets/d/1gooUY2Uh23d04bDt7Ph5gGQne4GB-LlApk5h1iO8aUA/edit#gid=0>`_.

Before 2023 there was little consistency, and they could be one of many adapter sets
available.
The full cornucopia of adapter sets available at the NGS facility is 
`here <https://ngs.vbcf.ac.at/forskalle3/account/adapters>`_, in a format that
could politely be called "a data-science nightmare".
Likely candidates are one of the eight "Nordborg Nextera INDEX" sets.
It may be best to ask within the lab in this case.

How can you determine which of the index sets was used?
The first place to look is `Forskalle <https://ngs.vbcf.ac.at/forskalle3/>`_.
Search (the bar is in the top right) for the sample number of your plate, then
look under "Adapter" under "Sample library options" in the upper left box, and
the index set information should be displayed.
For example, sample 144456 shows "Nordborg Nextera INDEX set 2..." in this box.

If in doubt, ask Viktoria to add the indices to the Nordborg-group NGS
submission sheet on google drive.

Work out which file is which
============================

From the previous section it is clear that if you know the adapter
sequence in the filename and which adapter set was used, you should be able to 
work out which file corresponds to which plate position. Doing this alignment is
fairly tedious, and there is no point in people replicating the task, so the
function ``align_fastq_with_plate_positions`` can done this automatically.
This function requires lists of input fastq files and a sample sheet giving the
details of each sample, and returns a copy of the sample sheet with the paths to
corresponding fastq files as additional columns

**It is important to note that this is not a foolproof procedure, and you need
to check the output.**

List the fastq files
--------------------

Using Python import the modules needed and make a list of filenames.
In this fictional example, imagine you have a folder of pairs of fastq files 
corresponding to forward and reverse read pairs for each sample. Files with mate
pair 1 and 2 are labelled ``R1`` and ``R2`` respectively. Globbing these files
returns a python list of files names.

.. code-block:: python

    from glob import glob
    import methlab as ml
    print("Using methlab version " + ml.__version__)
    # List of fastq files
    mate1=glob("path/to/bam_files/*_R1_*.fastq.gz")
    mate2=glob("path/to/bam_files/*_R2_*.fastq.gz")

Sequencing sample sheet 
-----------------------

To line these filenames up with biological samples we need a sample sheet, here
given as a CSV file and imported as a Pandas dataframe.

.. code-block:: python

    sequencing_sheet = pd.read_csv("tests/test_data/NGS_sample_sheet.csv")

Here is an example of how a sample sheet might look:

.. parsed-literal::
  
       plate      row  col    plantID  index_set      index1      index2
    0  2021-007   A    1  H1.2xH2.1-1          4  TTCTCGTGCA  ATACACAGAG
    1  2021-007   B    1  H1.2xH2.1-2          4  GCCTAACGTG  AGCTCTCAAG
    2  2021-007   C    1  H1.2xH2.1-3          4  CATTCACGCT  GTTGTACTCA
    3  2021-007   D    1  H1.2xH2.1-4          4  GCCATATAAC  ACACAATATC
    4  2021-007   D    1  H1.2xH2.1-4          4  TCGTGCATTC  GCGTTGGTAT
    5  2021-007   D    1  H1.2xH2.1-4          4  AACCAGCCAC  GCCACAGCAC

This file must include columns ``index1`` and ``index2`` that give the index
barcodes for each sample.
See the previous section for where to find these.
The other columns give information about each sample.
They can include any information, as long as ``index1`` and ``index2`` are present.
There may not be duplicate rows of ``index1`` and ``index2``, which in practical
terms mean that you should look at one sequencing plate at a time.

Be aware that sometimes index1 and index2 appear in reverse order in the file
names!

Merge filenames and sample sheet
--------------------------------

Pass the lists of file paths and the sample sheet to ``align_fastq_with_plate_positions``.

.. code-block:: python
  genotyping_sheet = ml.align_fastq_with_plate_positions(mate1, mate2, sequencing_sheet)    

The function adds two new columns called ``fastq1`` and ``fastq2`` giving the
paths from ``mate1`` and ``mate2`` for each sample it found a match for.

This function uses regular expressions to look for a nucleotide sequence inside
each filename and matches it to the indices in the sample sheet. This is a 
general, but **not very robust** solution. In particular, it will fail if there
are multiple samples with the same adpater sequences, such as when the same index
set was used for multiple plates. As such, it would be wise to **only use this on
one sequence plate at a time**. It also wouldn't hurt to look at the 
`source code <https://github.com/ellisztamas/methlab/blob/main/methlab/align_fastq_with_plate_positions.py>`
for this function and convince yourself of what it is doing.


Example
=======

Sequencing plate 2024-011 was 96 samples of F8 lines from the Intra-Swedish
collection. You can find the sequencing sample sheet, including plate positions,
sample names, barcode and other information
`here <https://docs.google.com/spreadsheets/d/1XOPm1GXGJ9zh1u43HGNsz9fKPbbBwL35tdAj-fhD0Vg/edit?usp=sharing>`.
I saved this information as a CSV file on the cluster as ``01_data/08_resequenced_F8s/sample_sheet.csv``.

The data were sequenced as NGS sample
`343561 <https://ngs.vbcf.ac.at/forskalle3/samples/343561>`.
After unzipping the data to ``/scratch-cbe/users/$(whoami)/04_resequencing/343561/``,
the first ten files look like this:

.. parsed-literal::

  343561_AACAAGTACAGGATATATCC_S87_I1_001.fastq.gz
  343561_AACAAGTACAGGATATATCC_S87_I2_001.fastq.gz
  343561_AACAAGTACAGGATATATCC_S87_R1_001.fastq.gz
  343561_AACAAGTACAGGATATATCC_S87_R2_001.fastq.gz
  343561_AACCGAGTTCGATCTCTGGA_S9_I1_001.fastq.gz
  343561_AACCGAGTTCGATCTCTGGA_S9_I2_001.fastq.gz
  343561_AACCGAGTTCGATCTCTGGA_S9_R1_001.fastq.gz
  343561_AACCGAGTTCGATCTCTGGA_S9_R2_001.fastq.gz
  343561_AACTTATCCTGTACTGGCGT_S78_I1_001.fastq.gz
  343561_AACTTATCCTGTACTGGCGT_S78_I2_001.fastq.gz

Below is an example Python script to link the sample names up to the file names.
If you haven't already, install the ``methlab`` and ``pandas`` modules for Python using 

.. code-block:: shell

  pip install methlab pandas

Here is the Python code:

.. code-block:: python

  import pandas as pd
  from glob import glob
  import os

  import methlab as ml
  print("Using methlab version " + ml.__version__)

  # === Input ===

  # Path to the directory containing raw fastq files
  path_to_fastq = "scratchdir/04_resequencing/343561/"
  if os.path.exists(path_to_fastq):
      print(f"Using fastq files from {path_to_fastq}")
  else:
      print(f"Error: could not find the path containing fastq files at {path_to_fastq}")
      print(f"Please check that you have unzipped the fastq files, and that you have the correct path.")
      print(f"The path is hard coded to a directory on the CLIP cluster, so you may need to change it.")
      print(f"Exiting...")
      exit()

  mate1 = glob( path_to_fastq + "/*_R1_*.fastq.gz")
  mate2 = glob( path_to_fastq + "/*_R2_*.fastq.gz")

  sequencing_sample_sheet = pd.read_csv(
      "01_data/08_resequenced_F8s/sample_sheet.csv"
  )

  genotyping_sample_sheet = ml.align_fastq_with_plate_positions(mate1, mate2, sequencing_sample_sheet)

``genotyping_sample_sheet`` should now contain all the columns in the original
sequencing sample sheet, plus additional columns with paths to the fastq files.