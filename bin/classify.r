library(lattice)

#----------------------------------------
# commandline args
#----------------------------------------

argv = commandArgs(TRUE)
plot_dir = argv[1]
input_gz_file = argv[2]
output_gz_file = argv[3]

#----------------------------------------
# global vars
#----------------------------------------

# exclude contigs less than this length
size.threshold = 2000

#----------------------------------------
# Load/filter/center/scale the data.
#
# Columns headers are: 
#
# 	1. contig (contig ID)
#	2. cov (mean read coverage)
#	3. gc (%gc content)
#	4. len (sequence length)
#----------------------------------------

data = read.table(gzfile(input_gz_file), header=TRUE)
data.subset = subset(data, len >= size.threshold)
data.subset.scaled = with(data.subset, data.frame(
		cov = scale(log(cov)),
		gc = scale(gc),
		len = scale(log(len))))

#----------------------------------------
# cluster the data with kmeans algorithm
#----------------------------------------

clusters = kmeans(data[c("cov","gc","len")], 2)
attach(clusters)

#----------------------------------------
# plot the clusters (for debugging purposes)
#----------------------------------------

pdf(paste(plot_dir, 'gc_vs_cov.pdf', sep='/'))
xyplot(cov ~ gc, data, group = cluster) 
invisible(dev.off())

pdf(paste(plot_dir, 'cov_vs_len.pdf', sep='/'))
xyplot(cov ~ len, data, group = cluster) 
invisible(dev.off())

pdf(paste(plot_dir, 'gc_vs_len.pdf', sep='/'))
xyplot(gc ~ len, data, group = cluster) 
invisible(dev.off())

#----------------------------------------
# determine which cluster contains the organelle
# contigs, by choosing the one with the higher
# mean coverage
#----------------------------------------

if (centers[1,"cov"] > centers[2,"cov"]) {
	organelle_cluster = 1
	genome_cluster = 2
} else {
	organelle_cluster = 2
	genome_cluster = 1
}

# change cluster index to recognizable category
cluster[cluster == organelle_cluster] = 'organelle'
cluster[cluster == genome_cluster] = 'genome'

#----------------------------------------
# Write out classification results. 
#
# Columns headers are input, but with
# an additional 'class' column: 
#
# 	1. contig (contig ID)
#	2. cov (mean read coverage)
#	3. gc (%gc content)
#	4. len (sequence length)
#	5. class ('organelle' or 'genome')
#----------------------------------------

data.out = cbind(data, cluster)
colnames(data.out) = c('contig', 'cov', 'gc', 'len', 'class')
file.out = gzfile(output_gz_file, 'w')
write.table(data.out, file.out, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
close(file.out)
