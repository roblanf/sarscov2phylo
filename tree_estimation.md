# What's the best method to get a tree of >>10K SARS-CoV-2 sequences?

Getting VERY large trees is difficult, but it's still useful to provide global context. And even if you're really only interested in the neighbourhood of a single sequence, it can be useful for figuring out what that sequences closest relatives are.

It's an open question (to me, at least) what the best method to get a good very large tree for SARS-CoV-2 sequences is. The point of this document is simply to list what I've tried and how I've compared them.

The basic methodology is as follows: 

1. Align the ~11K SARS-CoV-2 sequences from 2020-05-04-01 after trimming UTRs and removing a few bad sequences
1. Estimate a large tree in as many was as I practically can
2. Figure out what's best by comparing them in a likelihood framework
3. Also record execution times on 20 threads and 1 thread

I include command lines assuming you have a large (>>10K) alignment of GISAID SARS-CoV-2 sequences called global.fa

Here are the results. Methods are below.

## Results

### Which method is fastest? tl;dr rapidnj


| Tree inference method  | threads | time (s) |
|------------------------|---------|----------|
| quicktree default      | 1  	   | 1644	  |
| MASH->quicktree        | 20->1   | 2118 	  |
| IQ-TREE parsimony      | 20  	   | 1832	  |
| IQ-TREE parsimony      | 1  	   | 494	  |
| rapidnj_k2p      		 | 20  	   | 71	   	  |
| rapidnj_jc      		 | 20  	   | 72	   	  |
| rapidnj_k2p      		 | 1  	   | 186	  |
| rapidnj_jc      		 | 1  	   | 178	  |

So, IQ-TREE with parsimony on one thread, and `rapidnj` however you like to run it are a lot quicker than anything else. 

I've never seen anything as fast as `rapidnj`!

### Which method gives the best tree? tl;dr parsimony


Here I just took the tree from each of the five approaches I tried, ran it in IQ-TREE with a GTR and with/without re-estimating the branchlengths. The table shows the results parameters from each run. 

The result is fairly clear. The best method in both cases (with or without re-estimating branch lengths) is to estimate the tree with parsimony in IQ-TREE. This makes good sense - these are highly conserved genomes, and parsimony would be expected to do well. `rapidnj` also does very well, and the results are very close to those from the parsimony tree. The conclusion is that the best way to get a quick tree for large SARS-CoV-2 is to use parsimony. If you really need speed, then `rapidnj` is what you want. It will be worth keeping track of the relative speeds of parsimony and `rapidnj` as the dataset continues to increase in size. 

It's encouraging to see that the model parameters are very consistent, even though the tree topologies differ a lot. This means we can likely re-use these in later analyses to speed things up.

| method             | model | brlens       | lnL           | AICc        | delta_AICc  | tree_length | internal_brlen_pct | AC     | AG     | AT     | CG     | CT     | GT    |
|--------------------|-------|--------------|---------------|-------------|-------------|-------------|--------------------|--------|--------|--------|--------|--------|-------|
| quicktree default  | GTR   | fixed        | -217884.160   | 435784.331  | -175845.047 | 1.874       | 22.800             | 0.171  | 1.222  | 0.107  | 0.295  | 3.171  | 1.000 |
| MASH->quicktree    | GTR   | fixed        | -142614.8586  | 285245.7221 | -25306.438  | 0.4288      | 18.3429            | 0.2207 | 1.0489 | 0.1342 | 0.2471 | 3.0507 | 1.000 |
| IQ-TREE parsimony  | GTR   | fixed        | -129961.6397  | 259939.2843 | 0.000       | 0.288       | 20.769             | 0.217  | 0.981  | 0.139  | 0.196  | 3.002  | 1.000 |
| rapidnj_k2p        | GTR   | fixed        | -130224.1103  | 260464.2256 | -524.941    | 0.272       | 19.449             | 0.209  | 0.976  | 0.136  | 0.193  | 2.978  | 1.000 |
| rapidnj_jc         | GTR   | fixed        | -130231.2089  | 260478.4228 | -539.139    | 0.272       | 19.414             | 0.208  | 0.981  | 0.136  | 0.193  | 2.989  | 1.000 |
| quicktree default  | GTR   | re-estimated | -182757.4812  | 520138.5145 | -105602.523 | 0.451       | 17.688             | 0.171  | 1.220  | 0.107  | 0.296  | 3.148  | 1.000 |
| MASH->quicktree    | GTR   | re-estimated | -139291.3666  | 433206.2853 | -18670.294  | 0.318       | 19.842             | 0.220  | 1.050  | 0.134  | 0.247  | 3.053  | 1.000 |
| IQ-TREE parsimony  | GTR   | re-estimated | -129956.2197  | 414535.9915 | 0.000       | 0.288       | 20.723             | 0.217  | 0.981  | 0.139  | 0.196  | 3.002  | 1.000 |
| rapidnj_k2p        | GTR   | re-estimated | -129996.6322  | 414616.8164 | -80.825     | 0.288       | 21.016             | 0.211  | 0.987  | 0.138  | 0.193  | 3.010  | 1.000 |
| rapidnj_jc         | GTR   | re-estimated | -129996.0928  | 414615.7376 | -79.746     | 0.288       | 21.007             | 0.209  | 0.987  | 0.136  | 0.194  | 3.004  | 1.000 |


### Are the trees significantly different? tl;dr yes, don't use quicktree

Trees can differ, and on large alignments the likelihoods and AICc scores will look very different. But whether the differences are significant is another question. In short, we should only reject methods when the best tree is significantly better than the tree we get from those methods. To do this, we can use tree topology tests. 

The results of various tree topology tests are below. But the long and short of it is that, as one might predict from the data above, we can reject the `quicktree` trees, but the `rapidnj` and IQ-TREE trees are not significantly different from one another. That's good news, because it suggests that if we need speed it's OK to use `rapidnj` even though we might expect based on the likelihood that the IQ-TREE parsimony tree is marginally better.


|Tree | logL | deltaL  | bp-RELL | p-KH  | p-SH  | p-WKH  | p-WSH | c-ELW  | p-AU |
|-----|------|---------|---------|-------|-------|--------|-------|--------|------|
| quicktree default  |  -182866.3364 |   52910 |      0 - |     0 - |     0 - |     0 - |     0 -  |     0 -  | 8.59e-08 -  |
| MASH->quicktree    |  -139298.7327 |  9342.5 |      0 - |     0 - | 0.062 + |     0 - |     0 -  |     0 -  | 9.32e-108 - | 
| IQ-TREE parsimony  |  -129956.2194 |       0 |  0.522 + | 0.575 + |     1 + | 0.575 + | 0.859 +  | 0.522 +  |    0.571 +  |
| rapidnj_k2p        |  -129996.8756 |  40.656 |  0.227 + | 0.423 + | 0.839 + | 0.423 + | 0.875 +  | 0.227 +  |    0.421 +  |
| rapidnj_jc         |   -129996.426 |  40.207 |  0.251 + | 0.425 + | 0.832 + | 0.425 + | 0.865 +  | 0.251 +  |    0.443 +  |

deltaL  : logL difference from the maximal logl in the set.
bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990).
p-KH    : p-value of one sided Kishino-Hasegawa test (1989).
p-SH    : p-value of Shimodaira-Hasegawa test (2000).
p-WKH   : p-value of weighted KH test.
p-WSH   : p-value of weighted SH test.
c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002).
p-AU    : p-value of approximately unbiased (AU) test (Shimodaira, 2002).



### Which model is best

This is challenging, because estimating lots of parameters along with 22K branch lengths is going to be really hard. Nevertheless, it's worth it because if we can record the model parameters themselves, we can re-use these in later analyses since we don't expect them to change much as the trees change and the alignments grow (e.g. see above).

*results will be updated soon*


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

## Methods for comparing the trees

We have five trees to compare. There are lots of ways to compare trees, but a sensible one is to put them all in the same likelihood framework and just ask:

1. Which has the best likelihood (all assuming the same model)?
2. Can any reject each other?

Our five trees are:

* `quicktree_default.tree`
* `quicktree_mash.tree`
* `parsimony_iqtree.treefile`
* `rapidjn_from_aln_k2p.tree`
* `rapidjn_from_aln_jc.tree`

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

```

## Methods for topology tests
```
# first check that all trees end with a newline
sed -i -e '$a\' quicktree_default_nolb.tree
sed -i -e '$a\' quicktree_mash_nolb.tree
sed -i -e '$a\' parsimony_iqtree.treefile
sed -i -e '$a\' rapidnj_from_aln_k2p_noquotes.tree
sed -i -e '$a\' rapidnj_from_aln_jc_noquotes.tree

# now join our five trees of interest into one file
cat quicktree_default_nolb.tree quicktree_mash_nolb.tree parsimony_iqtree.treefile rapidnj_from_aln_k2p_noquotes.tree rapidnj_from_aln_jc_noquotes.tree > five_trees.tree

# check it's got five lines
wc -l five_trees.tree

# tree topology test, we'll use the starting tree as the parsimony tree (it's the best we've got) with optimised branch lengths
iqtree -s global.fa -z five_trees.tree -m GTR -te parsimony_iqtree.treefileGTR-varbr.treefile -nt 10 -keep-ident -zb 1000 -zw -au
```



## Methods for comparing models

We do each model indpendently here, because that allows us to also measure the execution times independently.

Of note, on earlier iterations of this analysis on a much smaller dataset (~2K seuqences) suggested that the +I+R3 model was the best. I add it here even though it will be expensive to estimate. The reason the free rate plus invariant sites models are likely to be good here is fairly simple. First, the invariant sites parameter accounts for the fact that a lot of the sites are constant, so will have effective rates of zero (the free-rate model itself doesn't allow for zero rates). The 3 rate categories of the R3 model from previous analyses included two low but similar rates, and one very very high rate. The super high rate most likely fits the data from the very very fast evolving sites (or are they common sequencing errors?) identified in this post from the Goldman lab: http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473. This is very useful, since if the inferential method can allow these sites to have very fast rates, they won't unduly harm the inference of the tree itself, and could on the other hand contribute to the inference of shallower parts of the tree, assuming that the sites are genuinely of biological origin, not due to technical errors. Of course, if one masks those sites as suggested in that post, then an R1 or R2 model may be a better fit to the data than the R3 model. But this analysis was initiated before that post was written, so the alignments I use here are not masked. Indeed, the only trimming for these alignments was to remove the UTRs.

```
tree=parsimony_iqtree.treefile
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I -nt 12 -me 0.05 -pre $tree'GTRI-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+G -nt 12 -me 0.05 -pre $tree'GTRG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+G -nt 12 -me 0.05 -pre $tree'GTRIG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 12 -me 0.05 -pre $tree'GTRIG-varbr' -redo


iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY -nt 12 -me 0.05 -pre $tree'HKY-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I -nt 12 -me 0.05 -pre $tree'HKYI-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+G -nt 12 -me 0.05 -pre $tree'HKYG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I+G -nt 12 -me 0.05 -pre $tree'HKYIG-varbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m HKY+I+R3 -nt 12 -me 0.05 -pre $tree'HKYIG-varbr' -redo


```





















