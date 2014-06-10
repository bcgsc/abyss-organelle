#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

# make 'sort' and 'join' use the same sort order
sort=LC_COLLATE=C sort

# Sub directories for intermediate files (to reduce clutter)
tmp=$(name)-data
plotdir=$(name)-plots

# number of threads
j?=1

# If "sortedsam" is provided, that SAM file is for the
# read-to-contig alignments, rather than generating them from
# scratch.

ifdef sortedsam
	sam_basename:=$(notdir $(shell basename -s .sam -s .sam.gz $(sortedsam)))
	bam:=$(sam_basename).bam
	sortedbam:=$(sam_basename).sorted.bam
endif
ifndef sortedbam
	bam:=$(name).bam
	sortedbam:=$(name).sorted.bam
endif

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check clean_aligns
default: args_check $(name).fa.gz

#------------------------------------------------------------
# argument checking
#------------------------------------------------------------

args_check:
ifndef ref
	$(error missing required arg 'ref')
endif
ifndef name
	$(error missing required arg 'name')
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

# make subdirs to store intermediate files (to reduce clutter)

$(tmp):
	mkdir -p $@

$(plotdir):
	mkdir -p $@

# index the reference

$(ref).bwt: $(ref)
	bwa index $(ref)

# compute per-contig coverage histograms

ifdef sortedsam
$(tmp)/$(name).genomecov.txt.gz: $(sortedsam) $(ref).bwt | $(tmp)
	smartcat $(sortedsam) | \
		samtools view -bSo - - | \
		bedtools genomecov -ibam - | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@
else
$(tmp)/$(name).genomecov.txt.gz: $(readfiles) $(ref).bwt | $(tmp)
	abyss-tofastq $(readfiles) | \
		bwa mem $(ref) - | \
		samtools view -bSo - - | \
		samtools sort -o - temp | \
		bedtools genomecov -ibam - | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@
endif

# compute per-contig coverage histograms to per-contig
# mean coverage

$(tmp)/$(name).cov.tab.gz: $(tmp)/$(name).genomecov.txt.gz
	zcat $< | cov-hist-to-mean | $(sort) | \
		(echo -e 'contig\tcov'; cat -) | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@

# compute %GC content per contig

$(tmp)/$(name).gc.tab.gz: $(ref)
	fastx2gc $(ref) | $(sort) | \
		(echo -e 'contig\tgc'; cat -) | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@

# compute contig lengths

$(tmp)/$(name).len.tab.gz: $(ref)
	bioawk -c fastx '{print $$name, length($$seq)}' $(ref) | $(sort) | \
		(echo -e 'contig\tlen'; cat -) | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@

# join contig ID, coverage, %GC, and length into one file

$(tmp)/$(name).cov_gc_len.tab.gz: \
		$(tmp)/$(name).cov.tab.gz \
		$(tmp)/$(name).gc.tab.gz \
		$(tmp)/$(name).len.tab.gz
	join -t $$'\t' <(zcat $(word 1, $^)) <(zcat $(word 2, $^)) | \
		join -t $$'\t' - <(zcat $(word 3, $^)) | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@

# classify contigs using kmeans clustering in R

$(tmp)/$(name).cov_gc_len_class.tab.gz: $(tmp)/$(name).cov_gc_len.tab.gz | $(plotdir)
	Rscript $(shell which classify.r) $(plotdir) $< $@

# build FASTA file of organelle contigs

$(name).fa.gz: $(tmp)/$(name).cov_gc_len_class.tab.gz $(ref)
	zcat $(word 1, $^) | \
		awk 'NR>1' | \
		join -t $$'\t' - <(bioawk -c fastx '{print $$name, $$seq}' $(word 2, $^) | $(sort)) | \
		$(sort) -n | \
		(echo -e 'contig\tcov\tgc\tlen\tclass\tseq'; cat -) | \
		bioawk -c header '$$class == "organelle" { print ">"$$contig; print $$seq }' | \
		gzip -c > $@.incomplete
	mv $@.incomplete $@
