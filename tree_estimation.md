# What's the best method to get a tree of >>10K SARS-CoV-2 sequences?

Getting VERY large trees is difficult, but it's still useful to provide global context. And even if you're really only interested in the neighbourhood of a single sequence, it can be useful for figuring out what that sequences closest relatives are.

It's an open question (to me, at least) what the best method to get a good very large tree for SARS-CoV-2 sequences is. The point of this document is simply to list what I've tried and how I've compared them.

The basic methodology is as follows: 

1. Estimate a large tree in as many was as I practically can
2. Figure out what's best by comparing them in a likelihood framework
3. Also record execution times on 20 threads and 1 thread

I include command lines assuming you have a large (>>10K) alignment of GISAID SARS-CoV-2 sequences called global.fa

Here are the results. Methods are below.

## Results

| Tree inference method  | threads | lnL (fixed brlen) | lnL (ML brlen) | time (s) |
|------------------------|---------|-------------------|----------------|----------|
| quicktree default      | 1  	   | 			 	   | 				| 1644	   |
| MASH->quicktree        | 20->1   |			       | 				| 2118 	   |
| IQ-TREE parsimony      | 20  	   |			       |				| 1832	   |
| IQ-TREE parsimony      | 1  	   |			       |				| 494	   |
| rapidnj_k2p      		 | 20  	   |			       |				| 71	   |
| rapidnj_jc      		 | 20  	   |			       |				| 72	   |
| rapidnj_k2p      		 | 1  	   |			       |				| 186	   |
| rapidnj_jc      		 | 1  	   |			       |				| 178	   |

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

### Comparing the trees

Now we have five trees to compare. There are lots of ways to compare trees, but a sensible one is to put them all in the same likelihood framework and just ask:

1. Which has the best likelihood (all assuming the same model)?
2. Can any reject each other?

Our five trees are:

* `quicktree_default.tree`
* `quicktree_mash.tree`
* `parsimony_iqtree.treefile`
* `rapidjn_from_aln_k2p.tree`
* `rapidjn_from_aln_jc.tree`

Let's compare them using IQ-TREE. A few initial analyses suggested here that more threads were a LOT faster, so I'll run all of these analyses on 40 threads. Note that we analyse every tree with and without optimising branch lengths. Also note that some of these anlayses require quite a bit of RAM (around 20GB).

Most published analyses I've seen of SARS-CoV-2 sequences use a GTR or a GTR+G model. I did a lot of model comparisons on the data a couple of weeks ago, and this suggested that GTR is the best model structure, but accounting for rate variation is best done with GTR+I+R3. (Next in line, in this order, were GTR+I+G, GTR+I, GTR+G, then GTR). The importance of the +I is not too surprising, most sites in the alignment do not have any segregating variation.

Preliminary analyses on this much larger alignment confirms that GTR+I+R3 is **by far** the best model here. So although it's slower to optimise by a fair bit, let's use it anyway because we're most interested in which of these trees is the best.


```
# quicktree trees are full of line breaks(!), so we fix that first
tr -d '\n' < quicktree_default.tree > quicktree_default_nolb.tree
tree=quicktree_default_nolb.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -fixbr -pre $tree'_fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -pre $tree'_varbr' -redo

tr -d '\n' < quicktree_mash.tree > quicktree_mash_nolb.tree
tree=quicktree_mash_nolb.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -fixbr -pre $tree'_fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -pre $tree'_varbr' -redo

# iqtree parsimony tree
tree=parsimony_iqtree.treefile
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -fixbr -pre $tree'_fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -pre $tree'_varbr' -redo

# rapidnj trees need to have quotes removed so IQ-TREE can read them
sed "s/'//g" rapidnj_from_aln_k2p.tree > rapidnj_from_aln_k2p_noquotes.tree
tree=rapidnj_from_aln_k2p_noquotes.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -fixbr -pre $tree'_fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -pre $tree'_varbr' -redo

sed "s/'//g" rapidnj_from_aln_jc.tree > rapidnj_from_aln_jc_noquotes.tree
tree=rapidnj_from_aln_jc_noquotes.tree
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -fixbr -pre $tree'_fixbr' -redo
iqtree -s global.fa -te $tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.05 -pre $tree'_varbr' -redo

```



iqtree -s global.fa -te rapidnj_from_aln_noquotes.tree -keep-ident -n 0 -m GTR+I+R3 -nt 40 -me 0.1 -pre egg_40thread -redo

























