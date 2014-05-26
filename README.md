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

# Usage

The main assembly script is `abyss-organelle.mk`. Usage of `abyss-organelle.mk` is the same as `abyss-pe` (see the [ABySS README.md](https://github.com/bcgsc/abyss#assembling-a-paired-end-library)), except that it requires an additional parameter:

* readfiles: a list of FASTA/FASTQ/SAM/BAM to align to the scaffolds in order to determine read coverage during the classification step

## Usage example

```
$ abyss-organelle.mk name=arabidopsis k=50 readfiles='reads1.fa.gz reads2.fa.gz' lib=readfiles
```

## Optional parameters

* classify: specify whether to do the classification at the contig stage or the scaffold stage. Possible values are `classify=scaffolds` or `classify=contigs` [`classify=scaffolds`]
* j: specify number of threads
