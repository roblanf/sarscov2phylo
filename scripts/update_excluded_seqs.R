#!/usr/bin/env Rscript

# usage: Rscript update_excluded_seqs.R old_iteration_file treeshrink_file
args = commandArgs(trailingOnly=TRUE)

treeshrink = args[1]
e = read.delim(treeshrink, header=F)
e = t(e)
e = e[complete.cases(e)]
f = data.frame(cbind(e, "on long branch identified by treeshrink q=0.05"))

write.table(f, file="excluded_sequences.tsv", sep="\t", col.names = F, row.names = F, quote=F)