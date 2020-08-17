# Moving to online phylogenetics

As the data keep growing in size, it makes less and less sense to start each analysis as if the previous analysis never happened. The current process acts like this, estimating every tree from scratch, even if the dataset only increases by a few percent in size.

The previous anlayses have shown that MP is an excellent method here, so in principle we could use this to place new sequences on a tree, and make some global and local improvements to that tree.

To do this, I wrote the script `global_tree_gisaid_start_tree.sh`. This script is very similar to the previous script in that it does the alignment from scratch. But this time, instead of using fasttree to make a NJ tree, which is then optimised with various ME and ML moves, we instead start with the tree from the previous iteration of the pipeline, then *add* the new sequences to it (there are various cleaning steps involved here, but the addition itself is done with MP in IQ-TREE), and finally optimise that tree with fasttree. It took me some time to get the steps of that pipeline together, but there are a number of questions that remain about the best way to optimise the starting tree.

The requirements are that it has to be *fast enough* to hopefully put out a new tree every day, and that we need some feeling of a guarantee that global and local issues in the tree could at least in principle have been ironed out. 

So, it seems sensible to include some global (i.e. SPRs with large radii) and local (lots of ML NNI) moves when optimising the tree. But I want to compare these (of course) against just estimating the tree from scratch as I usually would. Here are seven analyses which will allow me to answer a few questions that I set out below.

These analyses all use tree from 31-7-20 as the starting tree (where it's used) and GISAID data from 11-8-20 as the alignment. When using a starting tree, I take the following steps:

1. Remove sequences from the tree that are not in the alignment (this can happen because sequences are filtered)
2. Add to the tree sequences in the alignment that are not already in the tree (this uses IQ-TREE and MP, and took about 1hr to add about 3K sequences)

That's it. The total dataset here is about 52K sequences.


```
# 2 sprs, + ML NNIs
/usr/bin/time -o 1.mem.txt -v fasttree -nt -gamma -nni 0 -spr 2 -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 1.tree

# 0 sprs, + ML NNIs
/usr/bin/time -o 2.mem.txt -v fasttree -nt -gamma -nni 0 -spr 0 -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 2.tree

# 2 sprs, + ML NNIs, 5 rate categories
/usr/bin/time -o 3.mem.txt -v fasttree -nt -gamma -nni 0 -spr 2 -cat 5 -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 3.tree

# 2 sprs, + ML NNIs, 5 rate categories, no gamma
/usr/bin/time -o 4.mem.txt -v fasttree -nt -nni 0 -spr 2 -cat 5 -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 4.tree


# 2 sprs, + ML NNIs, 5 rate categories, no gamma, boot 100
/usr/bin/time -o 5.mem.txt -v fasttree -nt -nni 0 -spr 2 -cat 5 -boot 100 -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 5.tree

# 2 sprs, + ML NNIs, 5 rate categories, no gamma, no boot
/usr/bin/time -o 6.mem.txt -v fasttree -nt -nni 0 -spr 2 -cat 5 -nosupport -sprlength 200 -intree iqtree_seqsadded_mp.treefile $outputfasta > 6.tree


# no input tree, otherwise default except for increased SPR length
/usr/bin/time -o 7.mem.txt -v fasttree -nt -gamma -sprlength 200 $outputfasta > 7.tree
```

## Q1: What gives the best tree?

For this we just compare analysis 1 and analysis 7. Analysis 1 uses the previous iteration as a starting tree, where the new sequences were added by IQ-TREE.

Start tree lnL: -496120.668
No start tree lnL: -496332.143

Unsurprisingly, the analysis with the starting tree is better. It's not surprising because almost all the splits in this tree were *already* highly optimised. Adding a few thousand (about 6%) new sequences doesn't change many of the existing splits. But I add 2 global SPR moves (with ME) followed by fasttrees default number of ML NNIs (30 in this case). Of note, the analysis converged after just 6 rounds of ML NNIs, suggesting that the MP placements of the new sequences were very similar to the MP placements. That's encouraging. In addition, this is about 1 week of new data. So doing this every day should be very robust (and probably a little quicker). 

Of note the timings were also very different:

Start tree time: 10.17 hrs
No start tree time: 14.35 hrs

## Q2: How much time does each round of SPRs add?

Comparing analysis 1 (10.17 hrs) and analysis 2 (6.08 hrs), it's clear that each round of SPRs adds about 2 hrs to the analysis. That's about 30% of the analysis time without SPRs, so it's very significant. We should check if it's worth it...

## Q3: Are the SPR rounds worth it?

One round of global SPR moves is certainly worth it. If we don't do that we'll have zero rounds of SPR moves, and that means that any serious issues in the tree (introduced for whatever reason) will be very hard if not impossible to fix. NNI moves just don't cut it in this scenario. But it's still useful to know how much the 2 rounds of SPR moves improved the final log likelihood (also bearing in mind that the 2 SPR moves are done in an ME framework, so could in principle make things worse if ME optima conflict with ML optima)...

2 SPRs lnL: -496120.668
0 SPRs lnL: -496239.083

Surprisingly (to me at least) the two rounds of SPR moves really seemed to help. Seems worth investigating whether 1, 2, or 3 rounds of SPR moves is best in the next analysis.

## Q4: How much time can we save?

Analyses 2-6 are all about investigating ways of saving time. By comparing them to each other we can get an idea of how much time can be saved by changing various settings.

First, comparing 1 and 3 tells us if we save time or change the likelihood by moving from 20 to 5 rate categories:

cat 20 lnL: -496120.668
cat 5  lnL: -498323.359

Oof. That made the lnL a LOT worse. In addition, the cat 5 lnL analysis (number 3) still took 9.8hrs, so it is a very marginal time saving. Not worth it. 

Second, comparing 3 and 4 tells us if removing the -gamma flag speeds things up (it makes no difference to the tree, it just makes the log liklihoods comparable and scales the branch lengths a little, neither of which are all that important). Analysis 3 took 9.8 hrs, analysis 4 took 9.4 hrs, so it's about a 4% time saving to remove the -gamma flag. Not nothing, but not much.

Third, comparing 4 and 5 tells us if reducing the number of bootstraps from 1000 to 100 saves time. Analysis 4 took 9.4 hrs, and analysis 5 took 7.6 hrs. That's almost a 20% time saving. Seems worth it, particularly because 100 SH bootstraps is probably enough, and these numbers are very very imprecise anyway. 

Finally, comparing 5 and 6 tells us how much more we can save by removing bootstraps all together. Analysis 5 took 7.6 hrs, and analysis 6 took 7.2 hrs. So it's not a huge saving, and *some* measure of uncertainty in the splits seems worth having for a little bit of time. 

## Conclusions so far

The best settings are 20 rate categories, we can take or leave the `-gamma` flag, and we should use 100 bootstraps. 

## Some more investigation

The remaining question seems to be how many rounds of ME SPRs is best. So let's run a few additional analyses for that. Let's compare the likelihood and execution time for 0-4 rounds of SPRs, and for good measure let's try three spr lengths - 100, 200, 500.

```
/usr/bin/time -o 11.mem.txt -v fasttree -nt -gamma -nni 0 -spr 0 -sprlength 100 -log 11.log -intree iqtree_seqsadded_mp.treefile global.fa > 11.tree
/usr/bin/time -o 12.mem.txt -v fasttree -nt -gamma -nni 0 -spr 0 -sprlength 200 -log 12.log -intree iqtree_seqsadded_mp.treefile global.fa > 12.tree
/usr/bin/time -o 13.mem.txt -v fasttree -nt -gamma -nni 0 -spr 0 -sprlength 500 -log 13.log -intree iqtree_seqsadded_mp.treefile global.fa > 13.tree
/usr/bin/time -o 21.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 100 -log 21.log -intree iqtree_seqsadded_mp.treefile global.fa > 21.tree
/usr/bin/time -o 22.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 200 -log 22.log -intree iqtree_seqsadded_mp.treefile global.fa > 22.tree
/usr/bin/time -o 23.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 500 -log 23.log -intree iqtree_seqsadded_mp.treefile global.fa > 23.tree
/usr/bin/time -o 31.mem.txt -v fasttree -nt -gamma -nni 0 -spr 2 -sprlength 100 -log 31.log -intree iqtree_seqsadded_mp.treefile global.fa > 31.tree
/usr/bin/time -o 32.mem.txt -v fasttree -nt -gamma -nni 0 -spr 2 -sprlength 200 -log 32.log -intree iqtree_seqsadded_mp.treefile global.fa > 32.tree
/usr/bin/time -o 33.mem.txt -v fasttree -nt -gamma -nni 0 -spr 2 -sprlength 500 -log 33.log -intree iqtree_seqsadded_mp.treefile global.fa > 33.tree
/usr/bin/time -o 41.mem.txt -v fasttree -nt -gamma -nni 0 -spr 3 -sprlength 100 -log 41.log -intree iqtree_seqsadded_mp.treefile global.fa > 41.tree
/usr/bin/time -o 42.mem.txt -v fasttree -nt -gamma -nni 0 -spr 3 -sprlength 200 -log 42.log -intree iqtree_seqsadded_mp.treefile global.fa > 42.tree
/usr/bin/time -o 43.mem.txt -v fasttree -nt -gamma -nni 0 -spr 3 -sprlength 500 -log 43.log -intree iqtree_seqsadded_mp.treefile global.fa > 43.tree
/usr/bin/time -o 51.mem.txt -v fasttree -nt -gamma -nni 0 -spr 4 -sprlength 100 -log 51.log -intree iqtree_seqsadded_mp.treefile global.fa > 51.tree
/usr/bin/time -o 52.mem.txt -v fasttree -nt -gamma -nni 0 -spr 4 -sprlength 200 -log 52.log -intree iqtree_seqsadded_mp.treefile global.fa > 52.tree
/usr/bin/time -o 53.mem.txt -v fasttree -nt -gamma -nni 0 -spr 4 -sprlength 500 -log 53.log -intree iqtree_seqsadded_mp.treefile global.fa > 53.tree
```

I'll save those in a text file called spr.txt and run them like this:

```
parallel -j 15 -k --bar :::: spr.txt
```

| ID | SPRmoves | SPRlen | time     | lnL         | delta hrs | delta lnL |
|----|----------|--------|----------|-------------|-----------|-----------|
| 11 | 0        | 100    | 25613.75 | -496239.08  | 0.02      | 0.00      |
| 12 | 0        | 200    | 25551.77 | -496239.08  | 0.01      | 0.00      |
| 13 | 0        | 500    | 25529.44 | -496239.08  | 0.00      | 0.00      |
| 21 | 1        | 100    | 35276.59 | -496098.941 | 2.71      | 140.14    |
| 22 | 1        | 200    | 35207.22 | -496098.941 | 2.69      | 140.14    |
| 23 | 1        | 500    | 35239.38 | -496098.941 | 2.70      | 140.14    |
| 31 | 2        | 100    | 38495.2  | -496120.668 | 3.60      | 118.41    |
| 32 | 2        | 200    | 38688.24 | -496120.668 | 3.66      | 118.41    |
| 33 | 2        | 500    | 38522.25 | -496120.668 | 3.61      | 118.41    |
| 41 | 3        | 100    | 40822.96 | -496063.816 | 4.25      | 175.27    |
| 42 | 3        | 200    | 41290.1  | -496063.816 | 4.38      | 175.27    |
| 43 | 3        | 500    | 41456.04 | -496063.816 | 4.42      | 175.27    |
| 51 | 4        | 100    | 43605.31 | -495982.737 | 5.02      | 256.35    |
| 52 | 4        | 200    | 43502.2  | -495982.737 | 4.99      | 256.35    |
| 53 | 4        | 500    | 43689.79 | -495982.737 | 5.04      | 256.35    |

The best bang for our buck is to do just one round of SPRs (140 likelihood units per round, vs. e.g. 64 per round with 4 rounds). We get better likelihoods with more SPRs (in general) but the time penalty starts to get large. At 4 rounds of SPRs we've added 5 hrs to our analysis (meaning we're pushing our luck to get out an update every day, particularly as the dataset grows a lot in size).

Of note, I should have realised that an SPR radius of 100 encompasses the entire tree, so there is no difference in likelihood between any of the different SPR lengths. Still, I'll keep it at a very large number (1000) in future, to make sure that the SPR phase allows for truly global moves even as the tree grows a lot. 

So, to take a starting tree of ~50K SARS-CoV-2 sequences, add ~3K new sequences to it, and then optimise it in fasttree, the best commandline for that seems to be:

```
fasttree -nt -gamma -nni 0 -spr 1 -sprlength 1000 -boot 100 -intree iqtree_seqsadded_mp.treefile $outputfasta > 1.tree
```

Reasoning: I leave the gamma on there purely so that as things progress it's easy to compare the lnL of different runs; reduce bootstraps to 100, becuase the SH measures of support are potentially useful (though I'm not all that sure how to interpret them yet), but the additional precision of doing 1000 bootstaps is not justified, and use 1 spr move of length 1000 because that seems to be a good balance between likelihood and execution time. Though it's also tempting to consider having 4 SPR moves, because for a few additional hours it seems like we can get a higher likelihood. However, when this analysis runs every day (rather than, in this case, on 10 days worth of new sequences), I think we have a solution. 


There is one more thing I want to try, and that's whether we can make good improvements by running the 1 SPR analysis on the starting tree twice. The reason for this - it's often useful to interleave two algorithms (like SPR and NNI moves on trees) when optimising. So it's feasible that running an analysis with 1 SPR move followed by NNI moves twice, does a better job than any of the above, and could (let's see) give a nice balance of time and accuracy. 

To do this, we just want to compare this:

```
/usr/bin/time -o 91.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 1000 -boot 100 -intree iqtree_seqsadded_mp.treefile $outputfasta > 91.tree
```

to this:

```
/usr/bin/time -o 92.mem.txt -v fasttree -nt -nni 0 -spr 1 -sprlength 1000 -nosupport -intree iqtree_seqsadded_mp.treefile $outputfasta > 92.tree
/usr/bin/time -o 93.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 1000 -boot 100 -intree 92.tree $outputfasta > 93.tree
```

Note that in the second version, we can skip the gamma optimisation of branch lengths and the bootstraps, so it will less than double the execution time. And in the second step fo the second version, we use the tree from the first step as the starting tree.

Turns out this didn't really help much:

| Iterations | time     | lnL         |
|------------|----------|-------------|
| 1          | 22914.02 | -496098.941 |
| 2          | 39068.08 | -496038.076 |


I almost doubled the execution time, for a marginal gain in likelihood. 

## Final conclusion

This is the best solution so far. Iterated once a day with the new sequences from that day.

```
/usr/bin/time -o 91.mem.txt -v fasttree -nt -gamma -nni 0 -spr 1 -sprlength 1000 -boot 100 -intree iqtree_seqsadded_mp.treefile $outputfasta > 91.tree
```

P.S. In the end I decided to see what would happen if I just cranked it to do 10 SPRs, and the results were surprising. The lnL was -494195.750, which is a sizeable improvement. And the time was 11hrs. It seems worth keeping this for now, even if the 10 rounds of SPRs are not always worth it. 



