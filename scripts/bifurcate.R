#!/usr/bin/env Rscript

# input a tree with multifurcations
# output a tree with those multifurcations resolved

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ape))

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-i", "--input"), 
                     type = "character",
                     dest = 'input.tree',
                     help="Input tree"
                     )

parser <- add_option(parser, 
                     opt_str = c("-o", "--output"), 
                     type = "character",
                     dest = 'output.tree',
                     default=NA,
                     help="Output tree"
                     )

opt = parse_args(parser)


# bifurcate!
intree = read.tree(opt$input.tree)
outtree = multi2di(intree)
write.tree(outtree, file=opt$output.tree)