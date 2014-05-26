#!/usr/bin/make -Rrf

SHELL:=/bin/bash -o pipefail
PATH:=/home/benv/arch/genesis/bwa-0.7.4/bin:$(PATH)

# number of threads
j?=1

ifndef queryfiles
error::
	$(error missing parameter 'queryfiles')
endif
ifndef target
error::
	$(error missing parameter 'target')
endif
ifndef name
error::
	$(error missing parameter 'name')
endif

default: $(name).sorted.bam.bai

$(target).bwt: $(target)
	bwa index $(target)

$(name).sam.gz: $(target).bwt $(queryfiles)
	abyss-tofastq $(queryfiles) | bwa mem -t$j $(bwa_opt) $(target) - | gzip -c > $@

$(name).bam: $(name).sam.gz
	zcat $< | samtools view -bSo $@ -

$(name).sorted.bam: $(name).bam
	samtools sort $< $(name).sorted

$(name).sorted.bam.bai: $(name).sorted.bam
	samtools index $<
