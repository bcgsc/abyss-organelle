#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail
PATH:=/home/benv/arch/genesis/bedtools-2.17.0/bin:/home/benv/arch/genesis/bwa-0.7.4/bin:/home/benv/bin/myscripts:$(PATH)

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

ifdef sam
	bam:=$(notdir $(shell basename -s .sam -s .sam.gz $(sam))).sorted.bam
endif
ifndef bam
	bam:=$(name).sorted.bam
endif

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check
default: args_check $(name).genomecov.hist

#------------------------------------------------------------
# argument checking
#------------------------------------------------------------

args_check:
ifndef sam
ifndef ref
	$(error missing required arg 'ref')
endif
ifndef readfiles
	$(error missing required arg 'readfiles')
endif
endif
ifndef name
	$(error missing required arg 'name')
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

$(name).cov-per-contig.tab: $(name).genomecov.hist
	cov-hist-to-mean $< > $@

$(name).genomecov.hist: $(bam).bai
	bedtools genomecov -ibam $(bam) | egrep -v '^genome' > $@

ifdef sam
$(bam).bai: $(bam)
	samtools index $(bam)
$(bam): $(sam)
	smartcat $< | samtools view -bSo $(bam) -
else
$(bam).bai: $(ref) $(readfiles)
	bwa-mem.mk queryfiles='$(readfiles)' target='$(ref)' name='$(name)'
endif
