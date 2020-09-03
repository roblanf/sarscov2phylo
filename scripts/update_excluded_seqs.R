#!/usr/bin/env Rscript

# usage: Rscript update_excluded_seqs.R old_iteration_file treeshrink_file
args = commandArgs(trailingOnly=TRUE)

old_file = args[1]
old = read.delim(old_file, header=F)

treeshrink = args[2]
e = read.delim(treeshrink, header=F)
e = t(e)
e = e[complete.cases(e)]
f = data.frame(cbind(e, "on long branch identified by treeshrink q=0.05"))

names(f)=names(old)
new = rbind(old, f)

write.table(new, file="excluded_sequences.txv", sep="\t", col.names = F, row.names = F)