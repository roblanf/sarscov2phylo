# Phylogenetics for SARS-CoV-2

This repository contains scripts for reproducible phylogenetic analyses of SARS-CoV-2 data from GISAID.

For many analyses, running your own phylogenies won't be necessary - one can simply use the excellent NextStrain pipeline (particularly augur, for the phylogenetics).

We wanted to build something that did a little more than augur though, which means first solving a few problems in different ways, hence this repository. We hope that if the methods prove useful, we'll have time to build them into augur at a later date.

# Installation
You will need some version of conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

```
# 1. Clone the repo
git clone https://github.com/roblanf/sarscov2phylo``

# 2. change to the dir
cd sarscov2phylo/

# 3. download and install all the packages you need 
conda env create -f environment.yml

#4. acviate the conda environment with all the packages
conda activate sarscov2phylo``
```

# Quickstart

1. Get the latest GISAID data (full length and high coverage) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on GISAID.

2. Create the global alignment 

`bash gisad_to_global_alignment.sh -i gisaid.fasta -o global.fa -t 20 -k 100`

This runs a number of the auxilliary scripts, roughly doing the following things in the following order:

* fixes known issues with GISAID sequences
* excludes any sequences listed in `exclude_sequences.tsv`
* makes a reference alignment from the k (100 in the above example) most dissimilar sequences
* filters that alignment for bad sequences and sites

CURRENTLY NOT WORKING BUT SHOULD WORK SOON...
* profile aligns every sequence in the input fasta file to the reference alignment
* filters that alignment for bad columns and sequences to produce `global.fa`

On my server (not too flash) with 40 threads, this takes about 10 minutes to align a gisaid file with about 10K sequences in it. 

# Next steps under development

* build a global tree
* assign lineages with pangolin
* build lineage trees
* build high-quality local trees for a set of focal input genomes

# How you can help

* Have ideas for improving or adding to the methods: post on the issues and/or fork the repo and submit a pull request
* Know of SARS-CoV-2 genomes in GISAID that are problematic for any reason: post the sequence name and the reason on the issue tracker and I'll add them to the excluded_sequences.tsv file
 