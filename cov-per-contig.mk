#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail
PATH:=/home/benv/arch/genesis/bedtools-2.17.0/bin:/home/benv/arch/genesis/bwa-0.7.4/bin:/home/benv/bin/myscripts:$(PATH)

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

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
default: args_check $(name).cov_gc_len.tab.gz

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

$(name).cov_gc_len.tab.gz: $(name).genomecov.hist
	join -t $$'\t' \
		<(cov-hist-to-mean $< | sort) \
		<(bioawk -c fastx '{print $$name, gc($$seq), length($$seq)}' $(ref) | sort) | \
			sort -n | \
			gzip -c > $@

$(name).genomecov.hist: $(bam).bai
	bedtools genomecov -ibam $(bam) | egrep -v '^genome' > $@

ifdef sortedsam
$(bam).bai: $(bam)
	samtools index $(bam)
$(bam): $(sortedsam)
	smartcat $< | samtools view -bSo $(bam) -
else
$(bam).bai: $(ref) $(readfiles)
	$(if $(readfiles),,$(error missing required arg 'readfiles'))
	bwa-mem.mk queryfiles='$(readfiles)' target='$(ref)' name='$(name)'
endif
