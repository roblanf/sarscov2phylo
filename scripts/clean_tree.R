#!/usr/bin/env Rscript

# usage: Rscript clean_tree.R input_tree names_of_species_in_alignment.txt
# the names file should be one name per line

args = commandArgs(trailingOnly=TRUE)

library(ape)
library(readr)

print("Input tree")
print(args[1])
print("input names")
print(args[2])


t = read.tree(args[1])
a = readLines(args[2])

drop = t$tip.label[which(t$tip.label %in% a == FALSE)]

new = a[which(a %in% t$tip.label == FALSE)]

write_lines(new, "alignment_names_new.txt")

t2 = drop.tip(t, drop)

write.tree(t2, "input_tree_cleaned.tree")