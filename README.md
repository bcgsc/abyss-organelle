cov-per-contig.mk: do bwa-mem alignments and determine mean coverage per contig
cov-hist-to-mean: convert 'bedtools genomecov' coverage histogram output to "contig-name TAB mean-coverage"
fa-to-gc: print "contig-name TAB %gc-content" for each seq in FASTA file
