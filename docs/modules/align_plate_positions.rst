=================================
Align file names to 96-well plate
=================================

When we submit 96-well plates for sequencing to the NGS facility, we typically 
submit something like an Excel sheet giving rows/columns of the plate and the
biological sample that should be in each well. What we get back from the NGS
facility is a mass of bam or fastq files that look something like this::
    H3H7YDRXY_1#144456_ACTCGCTACGTCTAAT.bam

How is one meant to determine which sample in the original plate each file is
meant to correspond to?

What do the file names mean?
============================

The file name above contains information about the run which you can use to determine what has happened if you triangulate between the NGS master list and the NGS facility portal. Unfortunately, the naming system of files has changed through time, so they might not always look like the ones above. In this case at least you can see:

* ``H3H7YDRXY`` is the flow cell on which sequences were run.
* ``_1`` is some kind of subset of the data on that flow cell. For example, 1 and 2 here might indicate either end of paired-end data. Note that sometimes your data might be combined with someone else's data, so you might have 3 and 4, or some other complicated combination of data. It's best to ask Almudena or Viktoria if you aren't sure.
* ``144456`` is the facility sample number. You can use this to track down the facility's data on the sequencing run (in this case, for example: https://ngs.vbcf.ac.at/forskalle3/samples/144456). Confusingly, the facility also has a 'request number', which looks very similar.
* ``ACTCGCTACGTCTAAT`` gives the adapter index sequence for this sample.
This comprises two 8 nucleotide sequences that together give a unique identifier for the row/column position in a 96-well plate. There may be multiple combinations for separate plates so that these can be run on a ingle flow cell. For example, here is an example of the full set of Nextera Dual XT adapters, for up to four plates.
The full cornucopia of adapter sets available at the NGS facility is here, in a format that could politely be called "a data-science nightmare".