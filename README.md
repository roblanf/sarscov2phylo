# A global phylogeny of SARS-CoV-2

This repository provides a global phylogeny of SARS-CoV-2 data from GISAID, as well as reproducile methods with which to estimate it. It also provides code that can be used to estimate local phylogenies around your own sequences, please see below for more information.

We gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu™ Database on which this research is based. The list is detailed in [this tab-separated file](https://github.com/roblanf/sarscov2phylo/blob/master/gisaid_hcov-19_acknowledgement.tsv).

# Key files

The latest tree with the associated code to reproduce that tree can always be obtained from the latest release of this repository [here](https://github.com/roblanf/sarscov2phylo/releases/latest). Releases are named by the date (in Australian time) on which the sequences were downloaded from GISAID.

For convenience, you can also get the latest results via the following links:

* [Latest global phylogeny with Transfer Bootstrap Supports](https://github.com/roblanf/sarscov2phylo/blob/master/ft_TBE.tree)
* [Latest global phylogeny with Felsenstein Bootstrap Supports](https://github.com/roblanf/sarscov2phylo/blob/master/ft_FBP.tree)
* [Acknowledgements file for those that upload to GISAID](https://github.com/roblanf/sarscov2phylo/blob/master/gisaid_hcov-19_acknowledgement.tsv)
* [Latest version of the script to produce a global tree](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/global_tree_gisaid.sh)

Privacy rules around the alignments themselves mean that they cannot be released here. The alignments can be recreated by following the steps described below. If you are a GISAID member and would like a copy of the alignment, please email me and I'll share it with you.

Also, please note that the script to produce the global tree is continually modified as new data are released, so if you are looking at the master branch it may be the case that the latest code will not exactly reproduce the latest tree. If you would like to reproduce the latest tree, please head over to the [latest release](https://github.com/roblanf/sarscov2phylo/releases/latest), which contains the latest tree and the version of the code used to produce it.

![Part of global SARS-CoV-2 phylogeny](https://github.com/roblanf/sarscov2phylo/blob/master/tree_image.jpg)

# How do I look at these trees

The trees contain more than 30,000 sequences. If you would like to visualise them, I recommend using [Dendroscope](http://dendroscope.org/). This is very fast, has a very quick search that allows you to find sequences of interest, and allows you to switch easily between cladograms, phylogenies, etc. The image above is a screenshot from Dendroscope.


# Why are there two trees, and what are all the numbers?

The topology and branch lengths of the two trees are identical. In both cases, the topology is the best topology estimated by `fasttree` using the `-fastest` option, which [I found to perform the best out of a wide range of potential approaches](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md). The branch lengths represent substitutions per site. 

The only difference between the two trees are the support values represented for each branch. One tree (ft_FBP.tree) contains standard bootstrap values (also known as Felsenstein Bootstrap Proportions, or FBP). The other is built using the Transfer Bootstrap Expectation (or TBE). You can read more about these two values and how they differe [here](https://natureecoevocommunity.nature.com/users/87831-olivier-gascuel/posts/32426-renewing-felsenstein-s-phylogenetic-bootstrap-in-the-era-of-big-data) and [here](https://www.nature.com/articles/s41586-018-0043-0). 

The difference is fairly simple to describe. The FBP value represents the percentage of bootstrap trees in which a split appeared. The TBE value represents the proportion of taxa in a clade that are, on average, stable across bootstrap replicates. The TBE is always at least as high as the FBP. In trees of thousands of taxa (this tree has 10s of 1000s), FBP values can often be very close to zero simply because branches are very short and some taxa can move around a lot. The TBE is more informative in this case, because it measures how many taxa have to be moved to recover a branch in a tree, rather than just the presence/absence of that branch. 

# Why bother?

For most people this repository won't contain much of use. Instead, most people should use:

* [Nextstrain](https://nextstrain.org/ncov/global)
* [PANGOLIN](https://pangolin.cog-uk.io/)

This repository exists for two reasons:

1. To provide a global phylogeny of SARS-CoV-2 sequences that others might find useful
2. To develop and provide reproducible best-practice analyses for estimating that tree  

Figuring out the best approach for estimating a very large phylogeny like this is an open question. I have kept notes of methods I have tried and compared [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md). A big part of the reason for working on this is that my research group develops phylogenetic methods like [IQ-TREE](http://www.iqtree.org/), [PartitionFinder](http://www.robertlanfear.com/partitionfinder/), and [RWTY](https://cran.r-project.org/web/packages/rwty/vignettes/rwty.html). We are always working on improving these methods (indeed, some of these methods are used internally by Nextstrain and PANGOLIN), and as such we are working hard to make sure we understand these data as much as possible.

If you want to dig deeper into how we estimate the trees, and/or contribute by trying other methods or making suggestions, please read on. Feel free to use these trees for your work, but please cite this repository, and more importantly all the relevant tools that are used to make the trees. 

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

## Estimating a tree for just GISAID data (i.e. reproducing the trees in this repository)

To estimate a global phylogeny for the GISAID sequences, run [this script](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/global_tree_gisaid.sh) with the following commandline:

```
bash global_tree_gisaid.sh -i gisaid.fasta -o global.fa -t 20 -k 100
```

Here's what the code does:

1. Fixes known issues with GISAID sequences, include replacing spaces in sequences with N's, and illegal characters in names. (The EPI-ID is never changed, so if you are interested in a particular sequence, it's best to search for it with the EPI-ID)
2. Removes sequences that I or others have determined to be questionable (these are documented in [excluded_sequences.tsv](https://github.com/roblanf/sarscov2phylo/blob/master/excluded_sequences.tsv))
3. Trims the low-quality ends from every sequence in the input file (saved as `trimmed.fa`), following advice here: http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473.
4. Makes an alignment of the 100 most dissimilar (but high-quality, namely no more than 10 ambiguous bases, filtered with `seqmagick`) sequences in the input, then filters sites with >10% gaps from this alignment using `esl-alimask` (saved as `aln_k_filtered.fa`)
5. Creates a global alignment by aligning every sequence to `aln_k_filtered.fa` with MAFFT, parallelised with GNU parallel, and doing various format interconversions with the EASEL tools.
6. Filters sites from the alignment with >5 % gaps using `esl-alimask`.
7. Filters sequence from the alignment that are shorter than 28000 bp, or contain >1000 ambiguities using `esl-alimanip`
8. Estimates a global tree using `fasttree` with the `-fastest` setting ([why's this good?](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md))
9. Uses `goalign` and GNU `parallel` to make 100 bootstrap alignments, then make trees from these as in step 8.
10. Creates bootstrap support trees with standard bootstraps and the transfer bootstrap using `raxml-ng`. The latter measure is likely more appropriate for large SARS-CoV-2 trees (see above).
11. Removes sequences on very long branches from the tree using [TreeShrink](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4620-2) (these sequences are likely to be either of poor quality and/or poorly aligned, and are subsequently added to the list of excluded seuqences so they are not included in future iterations of the pipeline)
12. Roots the trees with [nw_utils](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2887050/) using sequence hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05, following advice [here](https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1)


## Estimating local trees for your own sequences

This is a slightly different problem from estimating a global tree. Instead of worrying about the global tree, we are now concerned with getting the best possible estimate of relationships around some focal sequence or sequences. To do this, the `global_tree_gisaid_plus.sh` does the following:

1. Follows steps 1-7 as above, with the exception that *none* of your focal sequences will be removed in the alignment filtering step. Some of their sites may be removed, but the sequence itself will be retained no matter its length or the number of N's it contains.
2. Estimates a global tree using `rapidnj`, which is very good at getting good trees very quickly ([why's this good?](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md))
3. Uses `goalign` and GNU `parallel` to estimate 100 bootstrap `rapidnj` trees.
4. For each of the focal sequences in your own data, the script then uses the `rapidnj` tree to select a sub-tree that contains at least 100 and at most 1000 of the most closely-related sequences, inbetween those limits, it will cut the tree at the first node that is at least `-d` (a user-supplied option) substitutions from the focal sequence. 
5. Re-estimates the local phylogeny around each focal sequence using `raxml-ng` with a GTR+I+G model of sequence evolution with 100 bootstraps. 
6. Creates bootstrap support trees with standard bootstraps and the transfer bootstrap using `raxml-ng`. The latter measure is likely more appropriate for large SARS-CoV-2 trees.
7. Attempts to make PDFs of each local tree with `Inkscape`, though you will need to install this yourself as the `conda` implementation of Inkscape seems to be dead right now.
 

# Key tools used

The scripts rely on a suite of excellent tools. Here are many of them:

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
* `nw_utils`