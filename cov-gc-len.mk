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

$(tmp):
	mkdir -p $(tmp)

$(name).cov_gc_len.tab.gz: $(tmp)/$(name).genomecov.hist
	join -t $$'\t' \
		<(cov-hist-to-mean $< | sort) \
		<(fa2gc $(ref) | sort) | \
		join -t $$'\t' - <(bioawk -c fastx '{print $$name, length($$seq)}' $(ref) | sort) | \
		sort -n | \
		(echo -e 'contig\tcov\tgc\tlen'; cat -) | \
		gzip -c > $@

$(tmp)/$(name).genomecov.hist: $(tmp)/$(bam).bai
	bedtools genomecov -ibam $(tmp)/$(bam) | egrep -v '^genome' > $@

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
