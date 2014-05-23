# Description

Utilities to assemble an organelle genome using ABySS.

# Dependencies

The following tools must be installed and available on your PATH:

* samtools
* bwa
* bioawk
* bedtools
* R

# Tools

* abyss-organelle.mk (main script): do a standard ABySS assembly and then extract organelle scaffolds
* classify.mk: split contigs file into 'organelle' and 'genome' contigs using k-means clustering on coverage, %GC content, and length
* classify.r: R script to do k-means clustering for classify.mk
* bwa-mem.mk: Makefile to do an alignment with bwa mem, sort, index, etc.
* cov-hist-to-mean: convert 'bedtools genomecov' output to mean coverage per contig
* fastx2gc: compute %GC content for each seq in a FASTA/FASTQ file
* smartcat: does zcat/bzcat/cat based on file type
