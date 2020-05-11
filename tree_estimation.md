# What's the best method to get a tree of >>10K SARS-CoV-2 sequences?

Getting VERY large trees is difficult, but it's still useful to provide global context. And even if you're really only interested in the neighbourhood of a single sequence, it can be useful for figuring out what that sequences closest relatives are.

It's an open question (to me, at least) what the best method to get a good very large tree for SARS-CoV-2 sequences is. The point of this document is simply to list what I've tried and how I've compared them.

The basic methodology is as follows: 

1. Align the ~11K SARS-CoV-2 sequences from 2020-05-04-01 after trimming UTRs and removing a few bad sequences
2. Estimate a large tree in as many was as I practically can
3. Figure out what's best by comparing them in a likelihood framework
4. Also record execution times on 20 threads and 1 thread
5. Take a look at what the best model is, and whether the tree can be improved via more thorough searching.

I include command lines assuming you have a large (>>10K) alignment of GISAID SARS-CoV-2 sequences called global.fa

Here are the results. Methods are below.

## Results

### Which method is fastest? 

tl;dr rapidnj

| Tree inference method  | threads | time (s) |
|------------------------|---------|----------|
| rapidnj_k2p      		 | 20  	   | 71	   	  |
| rapidnj_jc      		 | 20  	   | 72	   	  |
| rapidnj_jc      		 | 1  	   | 178	  |
| rapidnj_k2p      		 | 1  	   | 186	  |
| IQ-TREE parsimony      | 1  	   | 494	  |
| quicktree default      | 1  	   | 1644	  |
| IQ-TREE parsimony      | 20  	   | 1832	  |
| MASH->quicktree        | 20->1   | 2118 	  |
| fasttree fastest		 | 1 	   | 4061     |
| fasttree default 		 | 1	   | 4625	  |

These results speak for themselves. `rapidnj` is insanely quick. IQ-TREE random stepwise addition parsimony trees are also very quick. Everything else is a lot slower.

N.B. Here's a list of things I tried that didn't work:

* IQ-TREE with the -fast option, the BIONJ part of the inference is too slow (>1hr, and that's before we've started optimising it with ML)
* MPBoot: too slow for this data (i.e. doesn't add anything that IQ-TREE can't do, which is to estimate 1 tree in 8 minutes)
* raxml-ng: all attempts at using raxml-ng were too slow because even calcualting a single parsimony tree took >1hr

### Which method gives the best tree? 

tl;dr fasttree, but others are close

I calcualted the likelihood of each tree topology under a GTR model in IQ-TREE. I did it with and without re-estimating branch lengths, because that can help us figure out which method gives branch lengths that are closest to the branch lengths we get with ML approaches (which have some nicer properties than other methods for estimating branch lengths; results in the delta_lnL column).

Results in the table are ordered from best to worst by `delta_AICc`. `fasttree` is the best, and next best is a random stepwise addition tree in parsimony from IQ-TREE. `rapidnj` is very close beind parsimony. `quicktree` can be safely ignored as no use here. These results make sense - fasttree puts in a lot more time and effort to optimising the tree than the IQ-TREE parsimony tree or the `rapidnj` tree, both of which are one-shot trees without further optimisation.

| method             | lnL        | lnl_fixedbrlen | delta_lnL | AICc      | delta_AICc | tree_length |
|--------------------|------------|----------------|-----------|-----------|------------|-------------|
| fasttree_fastest   | -129223.52 | -129223.52     | 0.00      | 354507.62 | 0.00       | 0.27        |
| fasttree_default   | -129253.66 | -129253.67     | 0.00      | 354567.91 | -60.29     | 0.27        |
| IQ-TREE parsimony  | -129956.22 | -129961.64     | -5.42     | 414535.99 | -60028.37  | 0.29        |
| rapidnj_jc         | -129996.09 | -130231.21     | -235.12   | 414615.74 | -60108.11  | 0.29        |
| rapidnj_k2p        | -129996.63 | -130224.11     | -227.48   | 414616.82 | -60109.19  | 0.29        |
| MASH->quicktree    | -139291.37 | -142614.86     | -3323.49  | 433206.29 | -78698.66  | 0.32        |
| quicktree default  | -182757.48 | -217884.16     | -35126.68 | 520138.51 | -165630.89 | 0.45        |

### Are the trees significantly different? 

Trees can differ, and on large alignments the likelihoods and AICc scores will look very different. But whether the differences are significant is another question. In short, we should only reject methods when the best tree is significantly better than the tree we get from those methods. To do this, we can use tree topology tests. 

The results of various tree topology tests are below. The table shows p values for a range of tests, where a low p value means that tree is rejected under the test. All the tests make different assumptions, but as a rule they tend to be rather conservative (i.e. quick to reject trees).

Regardless, we can make some fairly clear conclusions. First, the fasttree trees and (probably) the parsimony tree are not significantly different from one another (the latter is rejected in many tests, but not by a huge amount, and this is a large tree and these are very conservative tests; on top of that the fasttree trees contain multifurcations which the parsimony tree doesn't). This is very encouraging, since parsimony was so much faster. The rapidnj trees are rejected by all but the SH test. Using a small set of trees like this is not really what the SH test is for, but a few things make me think that practically the rapidnj trees should probably not be ignored. First, these tests are conservative, so we should take their rejections with a pinch of salt. Second, the fasttree trees have a somewhat unfair advantage insofar as they contain multifurcations. All the other trees have minimum bounds on short branches, so their likelihoods will suffer considerably. 

All this is to say that I'd conclude that fasttree, parsimony, and rapidnj are all producing very good and very similar trees here. You might make different conclusions.


|Tree | logL | deltaL  | bp-RELL | p-KH  | p-SH  | p-WKH  | p-WSH | c-ELW  | p-AU |
|-----|------|---------|---------|-------|-------|--------|-------|--------|------|
| quicktree default  |  -182862.9989 |  53639  |     0 - |     0 - |     0 - |     0 - |     0 - |         0 - | 2.93e-68 -  | 
| MASH->quicktree    |  -139299.8604 |  10076  |     0 - |     0 - | 0.046 - |     0 - |     0 - |         0 - | 8.57e-06 -  | 
| IQ-TREE parsimony  |  -129956.6204 |  733.1  | 0.003 - | 0.002 - | 0.522 + | 0.002 - | 0.003 - |     0.003 - | 0.000679 -  | 
| rapidnj_k2p        |  -129996.7879 | 773.27  |     0 - |     0 - | 0.518 + |     0 - |     0 - |  3.1e-142 - | 7.68e-85 -  | 
| rapidnj_jc         |  -129996.2309 | 772.71  |     0 - |     0 - | 0.518 + |     0 - |     0 - | 7.14e-178 - | 1.07e-101 - |  
| fasttree_fastest   |  -129223.5178 |      0  | 0.625 + | 0.608 + |     1 + | 0.608 + | 0.969 + |     0.626 + |    0.655 +  |
| fasttree_default   |  -129253.7439 | 30.226  | 0.372 + | 0.392 + | 0.904 + | 0.392 + | 0.869 + |     0.371 + |    0.372 +  |


deltaL  : logL difference from the maximal logl in the set.
bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990).
p-KH    : p-value of one sided Kishino-Hasegawa test (1989).
p-SH    : p-value of Shimodaira-Hasegawa test (2000).
p-WKH   : p-value of weighted KH test.
p-WSH   : p-value of weighted SH test.
c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002).
p-AU    : p-value of approximately unbiased (AU) test (Shimodaira, 2002).


### Which model is best
tl;dr GTR+G is the best balance of execution time and model fit; GTR+I+R3 is the best.

This is challenging, because estimating lots of parameters along with 22K branch lengths is computationally very expensive. Nevertheless, it's worth it because if we can record the model parameters themselves, we can re-use these in later analyses since we don't expect them to change much as the trees change and the alignments grow (e.g. see above).

Models estimated in IQ-TREE using the parsimony tree.

| model    | AICc        | execution_time |
|----------|-------------|----------------|
| GTR+I+R3 | 408343.0549 | 25949.333      |
| GTR+I+G  | 409437.2883 | 47040.2        |
| GTR+G    | 410169.0518 | 1022.47        |
| GTR+I    | 410465.71   | 834.28         |
| GTR      | 414535.99   | 567.5          |


### Can tree search improve the fasttree tree

Both raxml-ng and IQ-TREE struggled because their starting trees took too long to estimate. We can potentially circumvent that by using the fasttree topology as a starting tree and asking if we can improve it with the more thorough search algorithms in IQ-TREE and raxml-ng.

Using IQ-TREE with either `-fast` or the default settings (GTR+G model, fixed parameters) did not improve the tree from the fasttree_fastest tree. The full IQ-TREE search took 30 hours on 4 threads. 

Using `raxml-ng` with default settings and a GTR+G model did improve the tree. 

* more results coming soon

### Viewing the tree

I asked on twitter how people like to view large trees, so I thought I may as well record what worked there too: https://twitter.com/RobLanfear/status/1258615595060178945

* Figtree - my go to. Works, but is not built for such huge datatsets so finding the taxa I'm interested in involves a lot of scrolling.

* Iroki (https://www.iroki.net/viewer): great promise, but just doesn't work. Lots of 'unresponsive' warnings from chrome. Took about half an hour just to *load* the tree. Didn't render it. Gave up.

* OneZoom (https://github.com/onezoom/OZtree): sounds like it would work well, but requires a lot of installation. Didn't try it. 

* 

## Methods for estimating large trees

This is all run on a single server, giving the programs that can be multithreaded 20 threads. Everything was also run on a single thread (except MASH, because it won't be quicker on 1 thread) to double check if lots of threads was actually helping, and because fast single-threaded runtimes will be useful for bootstrapping. 


### Method: quicktree default

Quicktree is fast, but it's single threaded. So we shouldn't expect miracles here. But let's see.

```
# first we convert to stockholm format
# I won't count this in the timing
esl-reformat stockholm global.fa > global.sto

start=`date +%s`

quicktree -in a -out t global.sto > quicktree_default.tree

end=`date +%s`

runtime=$((end-start))
echo $runtime

```

### Method: MASH distances and quicktree

This method could be tweaked in various ways to go faster. Also (of course) we don't actually need aligned sequences to get MASH distances. But alignments are necessary for downstream analyses, and I'm focussing here on figuring out the quickest and best way to go from alignmnet to tree.

Quicktree is single-threaded, but I use it here because it can take the output of `mash triangle` without modification. 

```
start=`date +%s`

mash triangle -k 21 -s 10000 -p 20 -i global.fa > global_mash_dists.phy
quicktree -in m -out t global_mash_dists.phy > quicktree_mash.tree

end=`date +%s`

runtime=$((end-start))
echo $runtime
```

### Method: parsimony with IQTREE

There are many reasons to think parsimony might perform well with these data. They are very highly conserved, and there will almost never (never say never) be more than one substitution on a single branch. Parsimony by random stepwise addition is also a very fast method.

The `-keep-ident` flag is important here. IQ-TREE usually removes identical seuqences and adds them back at the end. Nothing wrong with this, it's just a slow process that we can side-step. The other flags just stop IQ-TREE doing anything other than estimating a parsimony tree then quitting.

Here we use 1 thread because it turns out this is the fastest solution - too much crosstalk slows down the parsimony inference.

```
start=`date +%s`
iqtree -s global.fa -keep-ident -n 0 -m JC -fixbr -nt 1 -pre parsimony_iqtree
end=`date +%s`

runtime=$((end-start))
echo $runtime
```


### Method: rapidnj k2p

`rapidnj` is insanely fast. So fast one could certainly get bootstrap values on the tree, even with large trees, and this could really help in some cases. By default it uses the k2p model for distances.


```
start=`date +%s`
rapidnj global.fa -i fa -o t -c 20 -n > rapidnj_from_aln_k2p.tree
end=`date +%s`

runtime=$((end-start))
echo $runtime
```


### Method: rapidnj jc

By default rapidnj uses the k2p model. That's probably the best, but we can easily check like this.


```
start=`date +%s`
rapidnj global.fa -i fa -o t -c 20 -n -a jc > rapidnj_from_aln_jc.tree
end=`date +%s`

runtime=$((end-start))
echo $runtime
```


### Method: fasttree default

Fasttree might work here, let's see.

```
start=`date +%s`
fasttree -nosupport -nt global.fa > fasttree_default.tree
end=`date +%s`

runtime=$((end-start))
echo $runtime
```


### Method: fasttree fastest

Fasttree also has a `fastest` option.

```
start=`date +%s`
fasttree -nosupport -nt -fastest global.fa > fasttree_fastest.tree
end=`date +%s`

runtime=$((end-start))
echo $runtime
```



## Methods for comparing the trees

We have seven trees to compare. There are lots of ways to compare trees, but a sensible one is to put them all in the same likelihood framework and just ask:

1. Which has the best likelihood (all assuming the same model)?
2. Can any reject each other?

Our seven trees are:

* `quicktree_default.tree`
* `quicktree_mash.tree`
* `parsimony_iqtree.treefile`
* `rapidjn_from_aln_k2p.tree`
* `rapidjn_from_aln_jc.tree`
* `fasttree_default.tree`
* `fasttree_fastest.tree`

I compare them using IQ-TREE. A few initial analyses suggested here that more threads were a LOT faster, so I ran all of these analyses on 40 threads. Note that we analyse every tree with and without optimising branch lengths, and with various models. Also note that some of these anlayses require quite a bit of RAM (around 20GB).

Most published analyses I've seen of SARS-CoV-2 sequences use a GTR or a GTR+G model. I did a lot of model comparisons on the data a couple of weeks ago, and this suggested that GTR is the best model structure, but accounting for rate variation is best done with GTR+I+R3. The importance of the +I is not too surprising, most sites in the alignment do not have any segregating variation. Despite the difference in the model fit, early analyses suggested that these differences didn't significantly impact the topology, and GTR+I+R3 is a tough model to optimise on large trees, so we'll avoid that model for now and just work with the simpler variants of the GTR model since they're a lot faster to work with, and might allow us to refine the topology using ML in downstream analyses.

Since the trees will have different numbers of parameters depending on the model and whether we re-estimate branch lengths (of which there are ~22K), it's useful to compare them with an approach that takes account of the differing number of parameters, so I'll use the AICc. Note though that if the AICc of re-estimating branch lenghts is *worse* than that of the fixed branch length trees, this doesn't mean we shouldn't be estimating ML branch lenghts. When the branch lengths are fixed, we're still estimating them from the data - we just do it in the NJ stage. So although the AICc isn't counting them as free parameters in the ML analysis (correctly) they're still free parameters in the global sense. Because of that, I'll only compare things within the fixed vs. variable branch length categories, not between them. I'll also record model parameters to see how much they vary across analyses, simply because if we're able to fix model parameters in some downstream analyses that would be very useful for speeding things up.


```
# quicktree trees are full of line breaks(!), so we fix that first
tr -d '\n' < quicktree_default.tree > quicktree_default_nolb.tree
tree=quicktree_default_nolb.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

tr -d '\n' < quicktree_mash.tree > quicktree_mash_nolb.tree
tree=quicktree_mash_nolb.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

# iqtree parsimony tree
tree=parsimony_iqtree.treefile
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

# rapidnj trees need to have quotes removed so IQ-TREE can read them
sed "s/'//g" rapidnj_from_aln_k2p.tree > rapidnj_from_aln_k2p_noquotes.tree
tree=rapidnj_from_aln_k2p_noquotes.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

sed "s/'//g" rapidnj_from_aln_jc.tree > rapidnj_from_aln_jc_noquotes.tree
tree=rapidnj_from_aln_jc_noquotes.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

# fasttree trees
tree=fasttree_default.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo

tree=fasttree_fastest.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -fixbr -pre $tree'GTR-fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR -nt 40 -me 0.05 -pre $tree'GTR-varbr' -redo


```

## Methods for topology tests
```
# first check that all trees end with a newline
sed -i -e '$a\' quicktree_default_nolb.tree
sed -i -e '$a\' quicktree_mash_nolb.tree
sed -i -e '$a\' parsimony_iqtree.treefile
sed -i -e '$a\' rapidnj_from_aln_k2p_noquotes.tree
sed -i -e '$a\' rapidnj_from_aln_jc_noquotes.tree
sed -i -e '$a\' fasttree_fastest.tree
sed -i -e '$a\' fasttree_default.tree

# now join our five trees of interest into one file
cat quicktree_default_nolb.tree quicktree_mash_nolb.tree parsimony_iqtree.treefile rapidnj_from_aln_k2p_noquotes.tree rapidnj_from_aln_jc_noquotes.tree fasttree_fastest.tree fasttree_default.tree > all_trees.tree

# check it's got the right number of lines
wc -l all_trees.tree

# tree topology test, we'll use the starting tree as the best tree with optimised branch lengths
iqtree -s global.fa -z all_trees.tree -m GTR -te fasttree_fastest.treeGTR-varbr.treefile -nt 10 -keep-ident -zb 1000 -zw -au -pre topotests
```



## Methods for comparing models

We do each model indpendently here, because that allows us to also measure the execution times independently.

Of note, on earlier iterations of this analysis on a much smaller dataset (~2K seuqences) suggested that the +I+R3 model was the best. I add it here even though it will be expensive to estimate. The reason the free rate plus invariant sites models are likely to be good here is fairly simple. First, the invariant sites parameter accounts for the fact that a lot of the sites are constant, so will have effective rates of zero (the free-rate model itself doesn't allow for zero rates). The 3 rate categories of the R3 model from previous analyses included two low but similar rates, and one very very high rate. The super high rate most likely fits the data from the very very fast evolving sites (or are they common sequencing errors?) identified in this post from the Goldman lab: http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473. This is very useful, since if the inferential method can allow these sites to have very fast rates, they won't unduly harm the inference of the tree itself, and could on the other hand contribute to the inference of shallower parts of the tree, assuming that the sites are genuinely of biological origin, not due to technical errors. Of course, if one masks those sites as suggested in that post, then an R1 or R2 model may be a better fit to the data than the R3 model. But this analysis was initiated before that post was written, so the alignments I use here are not masked. Indeed, the only trimming for these alignments was to remove the UTRs.

```
tree=parsimony_iqtree.treefile
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I -nt 12 -me 0.05 -pre $tree'GTRI-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+G -nt 12 -me 0.05 -pre $tree'GTRG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+G -nt 12 -me 0.05 -pre $tree'GTRIG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 12 -me 0.05 -pre $tree'GTRIR3-varbr' -redo


iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY -nt 12 -me 0.05 -pre $tree'HKY-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I -nt 12 -me 0.05 -pre $tree'HKYI-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+G -nt 12 -me 0.05 -pre $tree'HKYG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I+G -nt 12 -me 0.05 -pre $tree'HKYIG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I+R3 -nt 12 -me 0.05 -pre $tree'HKYIR3-varbr' -redo


```


## Methods for updating the tree

Once we have a model and a tree, we can try updating the tree. But since the tree is so large, it will help to fix the model parameters first. We'll also use the best tree as our starting tree. I use four threads here as raxml-ng helpfully told me how to guess roughly the right number of threads for this data.

```
# try with fast
tree=fasttree_fastest.tree

# first we need to fix multifurcations in the fasttree tree
# we could probably just randomly resolve these, but this approach will use parsimony and a couple of rounds of optimisation with NNI
iqtree -s global.fa -g $tree -keep-ident -m 'GTR{0.1644,0.8062,0.1464,0.1433,2.7414}+G{0.2210}' -nt 4 -me 0.05 -pre fasttree_fastest_bifurcating -redo -fast

tree=fasttree_fastest_bifurcating.treefile
iqtree -s global.fa -t $tree -keep-ident -m 'GTR{0.1644,0.8062,0.1464,0.1433,2.7414}+G{0.2210}' -nt 4 -me 0.05 -pre improve-fast -redo -fast

# that didn't change the tree, but it's only two iterations.
iqtree -s global.fa -t $tree -keep-ident -m 'GTR{0.1644,0.8062,0.1464,0.1433,2.7414}+G{0.2210}' -nt 4 -me 0.05 -pre improve-slow -redo

# we'll also try raxml-ng
raxml-ng --msa global.fa --model GTR+G --threads 4 --tree $tree --prefix improve-raxml

```

## Methods for making a constraint tree

```

# first set all branches with length zero to have a bootstrap support of 0.999999
# note we also change their branchlength to be non-zero, simply to stop this sed
# command getting stuck in an infinite loop. In princple we could change it back later
# but why bother. 

# first make a tree with bootstraps
bash tree_mp.sh -i global.fa -t 50

# then set all minimum-length branches to have high TBE support
# minimum branch lenght in an IQ-TREE tree is 0.000001

tree=global.fa_mp_boot_TBE.tree
sed -e ':a' -e ' s/)[0-9]\+.[0-9]\+:0.000001/)0.999999:0.000000/g; ta' $tree > tree_a.tree

# now use nw_ed to collapse all branches with TBE<x
# we can choose x like this, noting that we are interested in the proportion of 
# nodes that we have fixed
nw_stats $tree # there are 21296 nodes in the tree
nw_ed $tree 'i & b<0.5' o | nw_stats - # if we fix all nodes >0.5 TBE we fix 16128 nodes! 

# so, fixing all minimum length branches and all those with TBE<0.5
# means we only have 21296-16128=5168 nodes to care about! That's nice.

nw_ed $tree 'i & b<0.5' o > constraint.tree

# now we optimise the tree by using this as a constraint tree and the other tree as a starting tree
# to see if this is any better at optimisation
# we'll use the MP tree as a starting tree, and the constraint tree to constrain the search
raxml-ng --msa global.fa --model GTR+G --threads 4 --tree-constraint constraint.tree --prefix raxml_constraint


# we'll compare directly by doing the same thing but without using the constraint tree
raxml-ng --msa global.fa --model GTR+G --threads 4 --tree $tree --prefix raxml_free


iqtree -s global.fa -t $tree -g constraint.tree -keep-ident -m GTR+G -nt 4 -pre iqtree_constraint

iqtree -s global.fa -t $tree -keep-ident -m GTR+G -nt 4 -pre iqtree_constraint

```
















