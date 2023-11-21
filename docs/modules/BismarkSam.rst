=====================
Disembowel a SAM file
=====================

`BismarkSam` is a class of tools to take apart a SAM file from the output of
Bismark when you need to mess with individual reads. In particular, this was
designed to get at the XM tag that indicates methylation status within a read.

Import a SAM file
=================

Read a sample SAM file into memory as a list of BismarkSam objects.
It is assumed that the SAM file came from a BAM file that had been aligned with
Bismark.

.. code-block:: python

    from methlab.BismarkSam import *
    path = "tests/test_data/chloroplast.sam"
    sam = read_SAM(path)

Return some basic information about the third read in the file: read ID, read
length, chromosome and read sequence.

.. code-block:: python

    sam[2].id
    sam[2].length
    sam[2].chr
    sam[2].seq

SAM tags
========

SAM files are just tab-delimited text files, where each entry says something
about the read.
We can access some of those.

Methylation of each cytosine
----------------------------

If the SAM file was made from a BAM that was aligned with Bismark it will
include the tag "XM:" which gives a string indicate metylation.
Pull out the XM tag with ``sam[2].xm_tag``.

The output looks something like this:

.. parsed-literal::

    'HHH.Z...Z..ZXH....H.....Z.Z..HH.H...............X.............hh........z.......z...................h....h............................'

- There is a character for each nucleotide in the read.
- ``.`` indicate a non-cytosine.
- ``h`` and ``H`` are CHH sites
- ``x`` and ``X`` are CHG sites
- ``z`` and ``Z`` are CG sites
- Capital letters indicate methylated cytosines; lower-case letters indicate unmethylated cytosines.

Strand
------

Different bisulphite protocols differ in how they process each of four strands.

 - OT: original top strand
 - CTOT: complementary to top strand
 - CTOB: complementary to bottom strand
 - OB: original botom strand

OT and OB strands are the original strands, and the complementary strands are
those synthesised during PCR. The latter are ignored in directional libraries,
but can be used in non-directional libraries.

If you run Bismark (from version 24.0 onward) with the option ``strandID`` it 
will include the additional tag ``YS:Z:`` in the output BAM/SAM file indicating 
strand (e.g. ``YS:Z:CTOT``). If present, ``BismarkSam`` includes this as the
attribute ``sam[2].strand``. If this tag is not present, it returns ``'NA'``.

Note that this information is difficult to retreive from the SAM bit flag, 
because Bismark changed their meanings! See this issue for discussion:
https://github.com/FelixKrueger/Bismark/issues/455

How many methylated sites?
==========================

* ``sam[2].count_mC()`` returns a list showing the total number of methylated and unmethylated cytosines
* ``sam[2].mC_per_read()`` returns a list showing the number of methylated reads, unmethylated reads, total read length, and whether or not cytosines are occur next to one another