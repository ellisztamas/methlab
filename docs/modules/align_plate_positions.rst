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

============================
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
  flow cell. For example, 
  `here <https://docs.google.com/spreadsheets/d/1gooUY2Uh23d04bDt7Ph5gGQne4GB-LlApk5h1iO8aUA/edit#gid=0>`_
  is an example of the full set of Nextera Dual XT adapters, for up to four plates.
  The full cornucopia of adapter sets available at the NGS facility is 
  `here <https://ngs.vbcf.ac.at/forskalle3/account/adaptors>`_, in a format that
  could politely be called "a data-science nightmare".

============================
Work out which file is which
============================

With ``methlab``
=======================

From the previous section it is clear that if you know the 16-nucleotide adapter
sequence in the filename and which adapter set was used, you should be able to 
work out which file corresponds to which plate position. Doing the alignment is
fairly tedious, and there is no point in people replicating the task, so the
function ``align_fastq_with_plate_positions`` can doe this automatically.

Using Python import the modules needed and make a list of filenames.
In this fictional example, imagine you have a folder of pairs of fastq files 
corresponding to forward and reverse read pairs for each sample.

.. code-block:: python

    import pandas as pd
    from glob import glob
    import os
    import methlab as ml
    print("Using methlab version " + ml.__version__)
    # List of fastq files
    input_files=glob("path/to/bam_files/*fastq")

Pass this list to the function ``align_fastq_with_plate_positions``.
This function looks for a nucleotide sequence inside each file name and matches
it to a dataframe giving row and column positions for that sequence.
It then returns a pandas data frame with a row for each sample, giving sample
name, and up to two file paths to the data files that correspond to it.
Single-end data will return a single path, paired-end data return two paths,
and more than two matches raise an exception. Note that this requires that the
names you give it come from a single plate, and will throw an error if there are
multiple samples matching the sample row/column combination.

Sample names are returned as the row/column name; here we specify an optional
prefix ``t4_p2_`` which adds the temperature and plate ID for these data, so
sample names will look something like ``t4_p2_G1``.

.. code-block:: python

    ml.align_fastq_with_plate_positions(
        input_files,
        adapter_indices = 'nordborg', # specify that we need the Nordborg group custom nextera indices
        prefix = "t4_p2_"
    )

Specifying the adpator index set
================================

You need to specify which adaptor index set was used.

Since we often use the same index sets over and over again, some commonly used sets are stored internally. These data were run with the 'Nordborg Nextera INDEX' set, which you can use by passing ``nordborg`` to adapter_indices. At the time of writing (July 2023) the lab is moving towards using a single set for everything, which is likely to be the 'Unique Nextera Dual XT' set, which you can use by passing ``dualxt`` to ``adapter_indices``. Note that I have don't have data using the dual XT adapters to test this with, so do some quality control on your data to double check it's working! 

You can specify a custom adaptor set by creating a CSV file giving each adaptor pair.
This file must contain at least the following columns, but may contain more:

- row : Row letter from A to H
- col : Column number from 1 to 12
- seq1 : Nucleotide sequence of index 1
- seq2 : Nucleotide sequence of index 2

If the index set contains combinations for more than one plate, include an 
additional column 'set'.

For example, here are the first 10 rows of the index file for the Unique Nextera
Dual XT index set:

.. parsed-literal::

    set,row,col,name1,seq1,name2,seq2
    1,A,1,7001,CGCTCAGTTC,5001,TCGTGGAGCG
    1,A,2,7002,TATCTGACCT,5002,CTACAAGATA
    1,A,3,7003,ATATGAGACG,5003,TATAGTAGCT
    1,A,4,7004,CTTATGGAAT,5004,TGCCTGGTGG
    1,A,5,7005,TAATCTCGTC,5005,ACATTATCCT
    1,A,6,7006,GCGCGATGTT,5006,GTCCACTTGT
    1,A,7,7007,AGAGCACTAG,5007,TGGAACAGTA
    1,A,8,7008,TGCCTTGATC,5008,CCTTGTTAAT
    1,A,9,7009,CTACTCAGTC,5009,GTTGATAGTG