# Phylogenetics for SARS-CoV-2

This repository contains scripts for reproducible phylogenetic analyses of SARS-CoV-2 data from GISAID.

For many analyses, running your own phylogenies won't be necessary - one can simply use the excellent NextStrain pipeline (particularly augur, for the phylogenetics).

I wanted to build something that did a little more than augur though, which means first solving a few problems in different ways, hence this repository. We hope that if the methods prove useful, we'll have time to build them into augur at a later date.

The point of this repository is to allow for reproducible phylogenetic analyses that more-or-less follow best practice. I had been somewhat frustrated by reading so many SARS-CoV-2 analyses that present sensible analyses but no code. This repo tries to follow best practice, and I appreciate all feedback which indicates how it could be improved. 

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
conda activate sarscov2phylo
```

# Quickstart

1. Get the latest GISAID data (full length and high coverage) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on GISAID.

2. Add your own sequences to the GISAID fasta file if you like, e.g. `cat gisaid.fa myseqs.fa -> allseqs.fa`

3. Run this: `bash sequences_to_tree.sh -i gisaid.fasta -o global.fa -t 20 -k 100`


# What does it do, exactly?

To see exactly what's happening, you can take a look at each of the scripts in the `/scripts` folder. They are all shell scripts and what they do is described at the top of each. 

Here's an overview of what the wrapper script `seqeucnes_to_tree.sh` is actually doing:

1. Fixes known issues with GISAID sequences (there are quite a few...)
2. Removes sequences that I or others have determined to be questionable 
3. Trims the UTRs from every sequence in the input file (saved as `trimmed.fa`)
4. Makes an alignment of the 100 most dissimilar (but high-quality) sequences in the input, then filters gappy sites from this (saved as `aln_k_filtered.fa`)
5. Creates a global alignment by aligning every sequence to `aln_k_filtered.fa` with MAFFT
6. Filters sequence from the alignment that are shorter than 29100 bp (that's ~1% shorter than the typical alignment), or contain >200 ambiguities (that's also ~1%)
7. Filters sites from the alignment with >1 % gaps.
8. Estimates a global tree with 100 bootstraps using `rapidnj`

On my server (which is not too flash) with 40 threads, the whole process takes a couple of hours. Most of that time is estimating the tree with bootstraps.

# Next steps under development

* build high-quality local trees for a set of focal input genomes

# How you can help

* Have ideas for improving or adding to the methods: post on the issues and/or fork the repo and submit a pull request
* Know of SARS-CoV-2 genomes in GISAID that are problematic for any reason: post the sequence name and the reason on the issue tracker and I'll add them to the excluded_sequences.tsv file
 
# Key tools used

The scripts rely on a suite of absolutely excellent tools. Here are many of them:

* `mafft`
* `rapidnj`
* Easel tools like `esl-alimanip`
* `iqtree`
* UCSC fasta tools like `faSomerecords`
* `seqmagick`
* GNU `parallel`
