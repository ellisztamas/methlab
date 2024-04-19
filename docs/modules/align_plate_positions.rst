"""""""""""""""""""""""""""""""""
Align file names to 96-well plate
"""""""""""""""""""""""""""""""""

When we submit 96-well plates for sequencing to the NGS facility, we typically 
submit something like an Excel sheet giving rows/columns of the plate and the
biological sample that should be in each well. What we get back from the NGS
facility is a mass of bam or fastq files that look something like this:

.. parsed-literal::

    H3H7YDRXY_1#144456_ACTCGCTACGTCTAAT.bam

How is one meant to determine which sample in the original plate each file is
meant to correspond to?

What do the file names mean?
============================

Unfortunately, the naming system of files has changed through time, so they might
not always look like the ones above. In this case at least you can see:

* ``H3H7YDRXY`` is the flow cell on which sequences were run.
* ``_1`` is some kind of subset of the data on that flow cell.
  For example, 1 and 2 here might indicate either end of paired-end data.
  Note that sometimes your data might be combined with someone else's data, so
  you might have 3 and 4, or some other complicated combination of data.
  It's best to ask Almudena or Viktoria if you aren't sure.
* ``144456`` is the facility sample number. You can use this to track down the
  facility's data on the sequencing run (in this case, for example: 
  https://ngs.vbcf.ac.at/forskalle3/samples/144456).
  Confusingly, the facility also has a 'request number', which looks very similar.
* ``ACTCGCTACGTCTAAT`` gives the **adapter index** sequence for this sample.
  This comprises two 8-or-more nucleotide sequences that together give a unique
  identifier for the row/column position in a 96-well plate. There may be
  multiple combinations for separate plates so that these can be run on a single
  flow cell.

Which adaptors do I have? 
=========================

For data created from 2023 or later you probably have the Nextera Dual XT
adapters, which are four sets of 96 unique pairs of adaptors.
You can view them 
`here <https://docs.google.com/spreadsheets/d/1gooUY2Uh23d04bDt7Ph5gGQne4GB-LlApk5h1iO8aUA/edit#gid=0>`_.

Before 2023 there was little consistency, and they could be one of many adaptor
available.
The full cornucopia of adapter sets available at the NGS facility is 
`here <https://ngs.vbcf.ac.at/forskalle3/account/adaptors>`_, in a format that
could politely be called "a data-science nightmare".
Likely candidates are one of the eight "Nordborg Nextera INDEX" sets.
It may be best to ask within the lab in this case.

How can you determine which of the indexs was used?
The first place to look is `Forskalle <https://ngs.vbcf.ac.at/forskalle3/>`_.
Search (the bar is in the top right) for the sample number of your plate, then
look under "Adaptor" under "Sample library options" in the upper left box, and
the index set information should be displayed.
For example, sample 144456 shows "Nordborg Nextera INDEX set 2..." in this box.

============================
Work out which file is which
============================

From the previous section it is clear that if you know the 16-nucleotide adapter
sequence in the filename and which adapter set was used, you should be able to 
work out which file corresponds to which plate position. Doing this alignment is
fairly tedious, and there is no point in people replicating the task, so the
function ``align_fastq_with_plate_positions`` can done this automatically.
This function requires lists of input fastq files and a sample sheet giving the
details of each sample, and returns a copy of the sample sheet with the paths to
corresponding fastq files as additional columns

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

Sample sheet 
------------

To line these filenames up with biological samples we need a sample sheet, here
given as a CSV file and imported as a Pandas dataframe.

.. code-block:: python

    sample_sheet = pd.read_csv("tests/test_data/NGS_sample_sheet.csv")

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

Merge filenames and sample sheet
--------------------------------

Pass the lists of file paths and the sample sheet to ``align_fastq_with_plate_positions``.

.. code-block:: python
  new_sheet = ml.align_fastq_with_plate_positions(mate1, mate2, sample_sheet)    

This function looks for a nucleotide sequence inside each filename and matches
it to the indices in the sample sheet.