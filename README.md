# A global phylogeny of SARS-CoV-2

This repository provides a regularly-updated global phylogeny of SARS-CoV-2 data from GISAID, as well as reproducile methods with which to estimate it. 

We gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu™ Database on which this research is based. The list is detailed in [this tab-separated file](https://github.com/roblanf/sarscov2phylo/blob/master/gisaid_hcov-19_acknowledgement.tsv).

# Key files

The latest tree with the associated code to reproduce that tree can always be obtained from the latest release of this repository [here](https://github.com/roblanf/sarscov2phylo/releases/latest). Releases are named by the date (in Australian time) on which the sequences were downloaded from GISAID.

For convenience, you can also get the latest results via the following links:

* [Latest global phylogeny with Transfer Bootstrap Supports](https://github.com/roblanf/sarscov2phylo/blob/master/ft_TBE.tree)
* [Latest global phylogeny with Felsenstein Bootstrap Supports](https://github.com/roblanf/sarscov2phylo/blob/master/ft_FBP.tree)
* [Latest global phylogeny with fasttree SH Supports](https://github.com/roblanf/sarscov2phylo/blob/master/ft_SH.tree)
* [Acknowledgements file for those that upload to GISAID](https://github.com/roblanf/sarscov2phylo/blob/master/gisaid_hcov-19_acknowledgement.tsv)
* [Latest version of the script to produce a global tree](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/global_tree_gisaid.sh)

Privacy rules around the alignments themselves mean that they cannot be released here. The alignments can be recreated by following the steps described below. If you are a GISAID member and would like a copy of the alignment for any specific tree in the releases, please email me and I'll share it with you.

Also, please note that the script to produce the global tree is continually modified as new data are released, so if you are looking at the master branch it may be the case that the latest code will not exactly reproduce the latest tree. If you would like to reproduce the latest tree, please head over to the [latest release](https://github.com/roblanf/sarscov2phylo/releases/latest), which contains the latest tree and the version of the code used to produce it.

![Part of global SARS-CoV-2 phylogeny](https://github.com/roblanf/sarscov2phylo/blob/master/tree_image.jpg)

# How do I look at these trees?

The trees contain more than 30,000 sequences. If you would like to visualise them, I recommend using [Dendroscope](http://dendroscope.org/). This is very fast, has a very quick search that allows you to find sequences of interest, and allows you to switch easily between cladograms, phylogenies, etc. The image above is a screenshot from Dendroscope.

# Why isn't the sequence I'm intersted in in the tree?

Sequences are filtered out for a few reasons:

1. Sequence is not marked as both 'complete' and 'high coverage' by GISAID
2. Shorter than 28,000 bp after alignment
3. Contains more than 1000 ambiguities after alignment
4. Identified as being on a long branch by TreeShrink
5. Identified as questionable in the [excluded_sequences.tsv](https://github.com/roblanf/sarscov2phylo/blob/master/excluded_sequences.tsv) file

If your sequence is in GISAID, and was submitted before the date noted in the [latest release of the repository](https://github.com/roblanf/sarscov2phylo/releases/latest), but it is not in the tree, then it was filtered for one of the above reasons.

# Why are there three trees, and what are all the numbers?

The topology and branch lengths of the three trees are identical. In all cases, the topology is the best topology estimated by `fasttree` with options tuned specifically for this dataset, see [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md) and [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation2.md). The branch lengths represent substitutions per site. You will see that if you multiply the branchlenghts by about 30,000 (which is roughly the length of the alignments) many of the branchlengths are close to integers. That's because there's very little variation in these sequences, meaning that many branches have some integer number of changes inferred on them.

The only difference between the three trees are the support values represented for each branch. One tree (ft_FBP.tree) contains standard bootstrap values (also known as Felsenstein Bootstrap Proportions, or FBP). The next has support values calculated using the Transfer Bootstrap Expectation (or TBE). You can read more about these two values and how they differe [here](https://natureecoevocommunity.nature.com/users/87831-olivier-gascuel/posts/32426-renewing-felsenstein-s-phylogenetic-bootstrap-in-the-era-of-big-data) and [here](https://www.nature.com/articles/s41586-018-0043-0). The third tree has support values calculated by the SH test in `fasttree`.

The difference between the TBE and FBP is fairly simple to describe. The FBP value represents the percentage of bootstrap trees in which a split appeared. The TBE value represents the proportion of taxa in a clade that are, on average, stable across bootstrap replicates. The TBE is always at least as high as the FBP. In trees of thousands of taxa (this tree has 10s of 1000s), FBP values can often be very close to zero simply because branches are very short and some taxa can move around a lot. The TBE is more informative in this case, because it measures how many taxa have to be moved to recover a branch in a tree, rather than just the presence/absence of that branch.

The measures are all designed to estimate a similar quantity - something related to the confidence one might have in a branch. But they take very different approaches to this. Here's a plot which shows what they look like on a tree of >30K SARS-CoV-2 sequences from 24/6/20. The plot shows clearly that FBPs tend to be the lowest, TBEs are higher, and SH values are even higher. FBPs and TBEs are well correlated, but neither correlates that well to the SH values from fasttree (there *is* a very significant correlation though, as you'd expect!). 

![A comparison of support measures](https://github.com/roblanf/sarscov2phylo/blob/master/support_image.jpeg)


# Why bother?

For most people this repository won't contain much of use. Instead, most people should use:

* [Nextstrain](https://nextstrain.org/ncov/global)
* [PANGOLIN](https://pangolin.cog-uk.io/)

This repository exists for two reasons:

1. To provide a global phylogeny of SARS-CoV-2 sequences that others might find useful
2. To develop and provide reproducible best-practice analyses for estimating that tree  

Figuring out the best approach for estimating a very large phylogeny like this is an open question. I have kept notes of methods I have tried and compared [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md), as well as sometimes in the issues of this repository. In short, all choices and decisions made at all stages of the analysis are evidnece based either from my own analyses or from published analyses of SARS-CoV-2 sequences. A big part of the reason for working on this is that my research group develops phylogenetic methods like [IQ-TREE](http://www.iqtree.org/), [PartitionFinder](http://www.robertlanfear.com/partitionfinder/), and [RWTY](https://cran.r-project.org/web/packages/rwty/vignettes/rwty.html). We are always working on improving these methods (indeed, some of these methods are used internally by Nextstrain and PANGOLIN), and as such we are working hard to make sure we understand these data as much as possible.

If you want to dig deeper into how we estimate the trees, and/or contribute by trying other methods or making suggestions, please read on. Feel free to use these trees for your work, but please cite this repository, and more importantly all the relevant tools that are used to make the trees. 

# Installation

You should be able to re-run the scripts I provide here, provided you have access to a reasonably large server (e.g. ~500GB RAM and about 50 cores). With those resources I can run the scripts with ~35 threads (as of early June 2020, each fasttree run takes about 10GB), and it will complete in around 24 hours.

To start, you will need some version of conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

Once you have that, getting all the pacakges you need is simple:

```
# 1. Clone the repo
git clone https://github.com/roblanf/sarscov2phylo

# 2. change to the dir
cd sarscov2phylo/

# 3. download and install all the packages you need 
conda env create -f environment.yml

# 4. acviate the conda environment with all the packages
conda activate sarscov2phylo
```

# Quickstart

1. Get the latest [GISAID data](https://www.gisaid.org/) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on [GISAID](https://www.gisaid.org/), though you will have to become a member to do this.

2. Run this: `bash global_tree_gisaid.sh -i gisaid.fasta -o global.fa -t 20`

You can get more information on what is happening from the scripts themselves, and below.

# What does it do, exactly?

Here's what the code does:

1. Fixes known issues with GISAID sequences. This includes: (i) ensuring all sequences start on a new line; (ii) excluding problematic sequences identified in a raft of different ways, listed and explained in [this file](https://github.com/roblanf/sarscov2phylo/blob/master/excluded_sequences.tsv); (iii) replaces spaces in sequence names with a "\_"; (iv) replaces brackets in a species name with "-"; (v) repalces spaces in sequences themselves (!!) with "N". Note that the EPI_ID in each sequence name still forms a unique identifier for each sequence, which can always be reconciled with data downloaded from GISAID. All of this is accomplished with the [clean_gisaid.sh](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/clean_gisaid.sh) script, and uses the following tools: `sed` and `faSomeRecords`. This creates the file ending with `_cleaned.fa`, where the start of the filename matches the name of the raw sequence data.

2. Creates a global sequence alignment. This is done by aligning every sequence to `NC_045512.2` from NCBI. This is achieved quickly by doing it one sequence at a time in parallel in the [global_profile_alignment.sh](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/global_profile_alignment.sh) script, then joining the individually-aligned sequences into a global alignment at the end. It uses the following tools: `mafft`, `faSplit`, `faSomeRecords`, `grep`, GNU `parallel`, and `find`. This creates the file `aln_global_unmasked.fa`.

3. Masks sites that are likely dominated by errors of various kinds from the global alignment. Specifically, all sites suggested in [this post](http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) are masked by replacing them with gaps, by downloading the latest (when the script is run) version of [this file](https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/subset_vcf/problematic_sites_sarsCov2.mask.vcf). Masking is performed with the [mask_alignment.sh](https://github.com/roblanf/sarscov2phylo/blob/master/scripts/mask_alignment.sh) script, which uses the following tools: `wget`, `faSplit`, `seqmagick`, `grep`, GNU `parallel`, `find`. This creates the file `aln_global_masked.fa`. Note that this also masks the first and last 30 informative base pairs of every sequence, as well as any 7bp window which contains at least 2 uninformative characters (gaps or N's).

4. Filters seuqences that are shorter than 28000 bp and/or have >1000 ambiguities. This is because both of these filters tend to be associated with lower quality sequences, and or those that tend to be missing sufficient data that it is hard to place them meaningfully in a phylogeny. This uses `esl-alimanip`, and creates the file `aln_global_filtered.fa`.

5. Removes sites that are >50% gaps, after converting N's to gaps. I convert N's to gaps because maximum likelihood phylogenetic analyses treat them as the same thing. I then remove all sites with >50% gaps (at the time of writing this has only ever been the sites masked in step 3) using `esl-alimask`. This creates the file `global.fa` (or whatever you specified as your output file name when running the script).

6. Calculates and prints to `alignments.log` the stats of all alignments, for simple sanity checking.

7. Estimates a global ML tree from `global.fa`. This is done by using `fasttree` with settings determined empirically to be the best, which are constantly updated see [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md) and [here](https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation2.md). The scripts then use `goalign` to create 100 bootstrap alignments followed by re-estatimating all the ML trees with `fasttree` as and the `-fastest` setting, using GNU `parallel` to manage parallelisaion. FBP and TBE values are calcualted with `gotree`, and SH values are calculted with `fasttree`. The resulting three trees are rooted with the seuqence 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' as suggested in [this preprint](https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1), using [`nw_reroot`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2887050/). This creates the files `global.fa_ft_TBE.tree`, `global.fa_ft_FBP.tree`, `global.fa_ft_SH.tree`.

8. Removes sequences on very long branches from the tree using [`TreeShrink`](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4620-2). These sequences are likely to be either of poor quality and/or poorly aligned, so rather unreliable to interpret in a phylogeny with such limited variation. They are subsequently added to the list of excluded seuqences so they are not included in future iterations of the pipeline. This creates the files `ft_TBE.tree` and `ft_FBP.tree`, and `ft_SH.tree`.

9. Roots the final two trees with [`nw_reroot`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2887050/) as in step 7. 


# One last thing

The trees and alignments produced by the scripts here are checked regularly by eye. No matter all the clever software tools, there is no replacement for the human eye when it comes to spotting anomalies. If I find an anomaly, I go back and change the scripts to deal with it automatically (e.g. by removing a seuqence or masking a site, changing alignment or tree building parameters, etc). Rest assured that each release contains the exact code that built the trees included in that release, and also that I would not publish a release with a known anomaly. Obviously I can't guarantee to find all of the problematic sequences or sites, but I try.
 