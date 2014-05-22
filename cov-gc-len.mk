#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail
PATH:=/home/benv/arch/genesis/bedtools-2.17.0/bin:/home/benv/arch/genesis/bwa-0.7.4/bin:/home/benv/bin/myscripts:$(PATH)

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

# Sub directories for intermediate files (to reduce clutter)
tmp=$(name)-tmp
plotdir=$(name)-plots

# If "sortedsam" is provided, that SAM file is for the
# read-to-contig alignments, rather than generating them from
# scratch.

ifdef sortedsam
	sam_basename:=$(notdir $(shell basename -s .sam -s .sam.gz $(sortedsam)))
	bam:=$(sam_basename).sorted.bam
endif
ifndef bam
	bam:=$(name).sorted.bam
endif

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check
default: args_check $(tmp)/$(name).cov_gc_len_class.tab.gz

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

# make subdir to store intermediate files (to reduce clutter)

$(tmp):
	mkdir -p $@

$(plotdir):
	mkdir -p $@

# build bam file for read-to-contig alignments

ifdef sortedsam
$(tmp)/$(bam).bai: $(tmp)/$(bam)
	samtools index $(tmp)/$(bam)
$(tmp)/$(bam): $(sortedsam) | $(tmp)
	smartcat $< | samtools view -bSo $(tmp)/$(bam) -
else
$(tmp)/$(bam).bai: $(ref) $(readfiles) | $(tmp)
	$(if $(readfiles),,$(error missing required arg 'readfiles'))
	bwa-mem.mk queryfiles='$(readfiles)' target='$(ref)' name='$(tmp)/$(name)'
endif

# use 'bedtools genomecov' to convert read-to-contig alignment to
# per-contig coverage

$(tmp)/$(name).genomecov.txt.gz: $(tmp)/$(bam).bai
	bedtools genomecov -ibam $(tmp)/$(bam) | \
		egrep -v '^genome' | \
		gzip -c > $@

# convert 'bedtools genomecov' output to mean coverage per contig

$(tmp)/$(name).cov.tab.gz: $(tmp)/$(name).genomecov.txt.gz
	zcat $< | cov-hist-to-mean | sort | \
		(echo -e 'contig\tcov'; cat -) | \
		gzip -c > $@

# compute %GC content per contig

$(tmp)/$(name).gc.tab.gz: $(ref)
	fa2gc $(ref) | sort | \
		(echo -e 'contig\tgc'; cat -) | \
		gzip -c > $@

# compute contig lengths

$(tmp)/$(name).len.tab.gz: $(ref)
	bioawk -c fastx '{print $$name, length($$seq)}' $(ref) | sort | \
		(echo -e 'contig\tlen'; cat -) | \
		gzip -c > $@

# join contig ID, coverage, %GC, and length into one file

$(tmp)/$(name).cov_gc_len.tab.gz: \
		$(tmp)/$(name).cov.tab.gz \
		$(tmp)/$(name).gc.tab.gz \
		$(tmp)/$(name).len.tab.gz
	join -t $$'\t' <(zcat $(word 1, $^)) <(zcat $(word 2, $^)) | \
		join -t $$'\t' - <(zcat $(word 3, $^)) | \
		sort -n | \
		gzip -c > $@

# classify contigs using kmeans clustering in R

$(tmp)/$(name).cov_gc_len_class.tab.gz: $(tmp)/$(name).cov_gc_len.tab.gz | $(plotdir)
	Rscript $(shell which classify.r) $(plotdir) $< $@
