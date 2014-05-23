#!/usr/bin/make -rRf

#------------------------------------------------------------
# environment
#------------------------------------------------------------

SHELL:=/bin/bash -o pipefail

#------------------------------------------------------------
# global vars
#------------------------------------------------------------

outdir=$(name)-organelle

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: args_check
default: args_check $(outdir)/$(name)-organelle.fa.gz

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

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

$(outdir):
	mkdir -p $@

$(name)-scaffolds.fa:
	abyss-pe scaffolds

$(outdir)/$(name)-organelle.fa.gz: $(name)-scaffolds.fa | $(outdir)
	classify.mk -C $(outdir) ref='../$<' readfiles='$(readfiles)' name='$(name)'
