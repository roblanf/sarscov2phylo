# A global phylogeny of SARS-CoV-2

This repository contains scripts that detail methods to estimate very large trees with bootstrap support (Felsenstein Bootstrap and Transfer Bootstrap) for SARS-CoV-2 data from GISAID.

Privacy rules around the alignments themselves mean that they cannot be released here. But the global phylogeny itself is available in [this tree file](https://github.com/roblanf/sarscov2phylo/blob/master/global_tree.phy). 

We gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu™ Database on which this research is based. The list is detailed in [this tab-separated file](https://github.com/roblanf/sarscov2phylo/blob/master/gisaid_hcov-19_acknowledgement.tsv).

# Why bother?

For most people this repository won't contain much of use. Instead, one should start with the amazing:

* [Nextstrain](https://nextstrain.org/ncov/global)
* [PANGOLIN](https://pangolin.cog-uk.io/)

This repository exists for two reasons:

1. To provide a global phylogeny of SARS-CoV-2 sequences that others might find useful
2. To develop and provide reproducible best-practice analyses for estimating that tree  

Figuring out the best approach for estimating a very large phylogeny like this is an open question. I have kept notes of methods I have tried and compared [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md). A big part of the reason for working on this is that my research group develops phylogenetic methods like [IQ-TREE](http://www.iqtree.org/), [PartitionFinder](http://www.robertlanfear.com/partitionfinder/), and [RWTY](https://cran.r-project.org/web/packages/rwty/vignettes/rwty.html). We are always working on improving these methods (indeed, some of these methods are used internally by the amazing Nextstrain and PANGOLIN services), so it's really important that my group really understands these data.

If you want to dig deeper into how we estimate the trees, and/or contribute by trying other methods or making suggestions, please read on. Feel free to use these trees for your work, but please cite this repository and more importantly all the relevant tools that are used to make the trees. 

# Installation

You should be able to re-run the scripts I provide here, provided you have access to a reasonably large server (e.g. ~500GB RAM and ~50 cores). 

To start, you will need some version of conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

Once you have that, getting all the pacakges you need is simple:

```
# 1. Clone the repo
git clone https://github.com/roblanf/sarscov2phylo``

# 2. change to the dir
cd sarscov2phylo/

# 3. download and install all the packages you need 
conda env create -f environment.yml

# 4. acviate the conda environment with all the packages
conda activate sarscov2phylo
```

# Quickstart

## Estimating a tree for just GISAID data

1. Get the latest [GISAID data](https://www.gisaid.org/) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on [GISAID](https://www.gisaid.org/), though you will have to become a member to do this.

2. Run this: `bash global_tree_gisaid.sh -i gisaid.fasta -o global.fa -t 20 -k 100`

You can get more information on what is happening from the scripts themselves, and below. 

## Estimating local trees for your own data

1. Get the latest [GISAID data](https://www.gisaid.org/) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on [GISAID](https://www.gisaid.org/), though you will have to become a member to do this.

2. Collect your own sequences into a FASTA file, e.g. `myseqs.fa`

3. Run this: `bash global_tree_gisaid_plus.sh -i gisaid.fasta -o global.fa -t 20 -k 100 -a myseqs.fa -d 2`

# What does it do, exactly?

Note that all trees are rooted with the WH4 sequence, as suggested here: 

## Estimating a tree for just GISAID data

To see exactly what's happening, you can take a look at each of the scripts in the `/scripts` folder. They are all shell scripts and what they do is described at the top of each. 

1. Fixes known issues with GISAID sequences (there are quite a few...)
2. Removes sequences that I or others have determined to be questionable (these are documented in [excluded_sequences.tsv](https://github.com/roblanf/sarscov2phylo/blob/master/excluded_sequences.tsv))
3. Trims the low-quality ends from every sequence in the input file (saved as `trimmed.fa`), following the suggestions in [this virological post](http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473)
4. Makes an alignment of the 100 most dissimilar (but high-quality) sequences in the input, then filters gappy sites from this (saved as `aln_k_filtered.fa`)
5. Creates a global alignment by aligning every sequence to `aln_k_filtered.fa` with MAFFT
6. Filters sequence from the alignment that are shorter than 29100 bp (that's ~1% shorter than the typical alignment), or contain >200 ambiguities (that's also ~1%)
7. Filters sites from the alignment with >1 % gaps.
8. Estimates a global tree using `fasttree` with the `-fastest` setting ([why's this good?](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md))
9. Uses `goalign` and GNU `parallel` to estimate 100 bootstrap trees as in step 7.
10. Creates bootstrap support trees with standard bootstraps and the transfer bootstrap using `raxml-ng`. The latter measure is likely more appropriate for large SARS-CoV-2 trees.


## Estimating local trees for your own sequences

This is a slightly different problem. Instead of worrying about the global tree, we are now concerned with getting the best possible estimate of relationships around some focal sequence or sequences. To do this, the `global_tree_gisaid_plus.sh` does the following:

1. Follows steps 1-7 as above, with the exception that *none* of your focal sequences will be removed in the alignment filtering step. Some of their sites may be removed, but the sequence itself will be retained no matter its length or the number of N's it contains.
2. Estimates a global tree using `rapidnj`, which is very good at getting quite good trees ([why's this good?](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md))
3. Uses `goalign` and GNU `parallel` to estimate 100 bootstrap `rapidnj` trees.
4. For each of the focal sequences in your own data, the script then uses the `rapidnj` tree to select a sub-tree that contains at least 100 and at most 1000 of the most closely-related sequences, inbetween those limits, it will cut the tree at the first node that is at least `-d` (a user-supplied option) substitutions from the focal sequence. 
5. Re-estimates the local phylogeny around each focal sequence using `raxml-ng` with a GTR+I+G model of sequence evolution with 100 bootstraps. 
6. Creates bootstrap support trees with standard bootstraps and the transfer bootstrap using `raxml-ng`. The latter measure is likely more appropriate for large SARS-CoV-2 trees.
7. Attempts to make PDFs of each local tree with `Inkscape`, though you will need to install this yourself as the `conda` implementation of Inkscape seems to be dead right now.
 

# Key tools used

The scripts rely on a suite of absolutely excellent tools. Here are many of them:

* `mafft`
* `raxml-ng`
* `iqtree`
* `fasttree`
* `rapidnj`
* Easel tools like `esl-alimanip`
* UCSC fasta tools like `faSomerecords`
* `seqmagick`
* GNU `parallel`
* `conda`
* `goalign`