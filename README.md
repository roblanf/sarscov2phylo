# Phylogenetics for SARS-CoV-2

This repository contains scripts for reproducible analyses of SARS-CoV-2 data using GISAID

For many analyses, running your own phylogenies won't be necessary - one can simply use the excellent NextStrain pipeline (particularly augur, for the phylogenetics).

We wanted to build something that did a little more than augur though, which means first solving a few problems in different ways, hence this repository. We hope that if the methods prove useful, we'll have time to build them into augur at a later date.


# Quickstart
Get the latest GISAID data (full length and high coverage) in a single fasta file, e.g. `gisaid.fa`

`sh anlaysis.sh -i gisaid.fa -o tree.phy -t 8 -k 100`

So far this should:

* create an alignment of the k=100 most divergent sequences

Next steps

* fix profile-alignment of the rest of the seuqences to this alignment 
* filter the global alignment for bad sites and sequences
* build a global tree