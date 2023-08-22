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