#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

outdir=$(name)-organelle
classify?=scaffolds

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check
default: args_check $(outdir)/$(name)-organelle-scaffolds.fa

#------------------------------------------------------------
# argument checking
#------------------------------------------------------------

args_check:
ifndef readfiles
	$(error missing required arg 'readfiles')
endif
ifndef name
	$(error missing required arg 'name')
endif
ifeq ($(shell echo $(classify) | egrep '^(contigs)|(scaffolds)$$'),)
	$(error please use either 'classify=contigs' or 'classify=scaffolds' (default))
endif

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

$(outdir):
	mkdir -p $@

ifeq ($(classify),scaffolds)

# classify=scaffolds

$(name)-scaffolds.fa:
	abyss-pe scaffolds

$(outdir)/$(name)-organelle-scaffolds.fa: $(name)-scaffolds.fa | $(outdir)
	classify.mk -C $(outdir) ref='../$<' readfiles='$(readfiles)' name='$(name)-organelle-scaffolds'
	gunzip $(outdir)/$(name)-organelle-scaffolds.fa.gz

else

# classify=contigs

$(name)-6.fa:
	abyss-pe contigs

$(outdir)/$(name)-organelle-6.fa: $(name)-6.fa | $(outdir)
	classify.mk -C $(outdir) ref='../$<' readfiles='$(readfiles)' name='$(name)-organelle-6'
	gunzip $(outdir)/$(name)-organelle-6.fa.gz

$(outdir)/$(name)-organelle-6.dot: $(outdir)/$(name)-organelle-6.fa
	abyss-overlap -v $< > $@

$(outdir)/$(name)-organelle-scaffolds.fa: \
		$(outdir)/$(name)-organelle-6.fa \
		$(outdir)/$(name)-organelle-6.dot
	abyss-pe -C $(outdir) name='$(name)-organelle' \
		-o $(name)-organelle-6.fa -o $(name)-organelle-6.dot \
		scaffolds

endif
