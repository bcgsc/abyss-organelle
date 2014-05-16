#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail
PATH:=/home/benv/arch/genesis/bedtools-2.17.0/bin:/home/benv/arch/genesis/bwa-0.7.4/bin:/home/benv/bin/myscripts:$(PATH)

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check
default: args_check $(name).genomecov.txt

#------------------------------------------------------------
# argument checking
#------------------------------------------------------------

args_check:
ifndef ref
	$(error missing required arg 'ref')
endif
ifndef readfiles
	$(error missing required arg 'readfiles')
endif
ifndef name
	$(error missing required arg 'name')
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

$(name).cov-per-contig.tab: $(name).cov-per-contig.hist
	cov-hist-to-mean $< > $@

$(name).cov-per-contig.hist: $(name).sorted.bam.bai $(readfiles)
	bedtools genomecov -ibam $(name).sorted.bam | egrep -v '^genome' > $@

$(name).sorted.bam.bai: $(ref) $(readfiles)
	bwa-mem.mk queryfiles='$(readfiles)' target='$(ref)' name='$(name)'
