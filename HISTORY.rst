=======
History
=======

v0.6.0 (21st November 2023)
---------------------------
Changed the package name to methlab.

v0.5.0 (13th November 2023)
---------------------------
Added ``methylation_state`` to quantify the evidence that regions are unmethylated, CG-only methylated or TE-like methylated.

v0.4.3 (25th October 2023)
--------------------------
Added support for checking methylation in windows for only a subset of chromosomes

v0.4.2 (15th September 2023)
----------------------------
- Committed docs for ``CytosineCoverageFile`` and ``BismarkSam``.
- Added support for pulling out strand tags from SAM files. This generalises to
    any custom SAM tag. 

v0.4.1 (12th September 2023)
----------------------------

- Fixed ``install_requires``.
- SAM files now return the SAM flag.


v0.4.0 (12th September 2023)
----------------------------

- Class to manipulate aligned SAM files from Bismark.  
- Included docs on how to use.
    

v0.3.0 (8th September 2023)
---------------------------

* Documentaion for align_fastq_with_plate_positions and CytosineCoverageFile.

v0.2.2 (25th July 2023)
-----------------------

Created a class to import and take apart cytosine coverage files from Bismark.

- `count_reads` counts up methylated and unmethylated reads, plus total number of cytosines in each sequence context.
- `methylation_over_features` takes series of genome coordinates (from an
    annotation file, for example), pulls out each part of the coverage file and
    calls `count_reads` on each.
* `conversion_rate` calculates mean methylation on each chromosome
* `methylation_in_windows` partitions the genome in to windows of fixed size and
    calls `count_reads` on each.

0.1.0 (2023-06-29)
------------------

- Initial commit 
- Function `align_fastq_with_plate_positions.py`` looks up adaptor index positions from bam file names
