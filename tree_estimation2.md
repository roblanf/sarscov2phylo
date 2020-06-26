# What fasttree settings give the best tree of >>40K SARS-CoV-2 sequences?

The data have grown a lot, from ~10K to ~40K sequences, since I last compared methods to estimate trees. Since then I've been succesfully using `fasttree`. However, using data on the 14th of June 2020, I and others noticed that one clade in the tree is placed in what is very likely to be the wrong place (it's clade A.3 according to the PANGOLIN notation, if you are interested, and it's noted in the release for that date, because it sticks out a like a sore thumb when you look at the tree).

So, it's time to revisit fasttree to see if I can run it better and avoid the misplacement of this clade. I'll use this as an opportunity to run through a whole range of `fasttree` settings on this alignment, because it's obviously challenging.

It's not possible to search through all parameter combinations, there are just too many. So my strategy is to do a greedy parameter search. I.e. I'll start with the settings I have been using, then try other settings a few at a time, then just keep the best settings as I go. Sometimes I might revisit old settings. Think of it as hill climbing guided by human intuition and some evidence from watching fasttree work a lot. The nice thing is this - as long as I keep the alignment the same (and I will) I can just keep choosing the settings that give me the best likelihood (or AICc score, if the number of parameters is different). In each iteration, I'll usually only vary things by one setting from the settings in the previous iteration. Finally, note that we care about runtime as well as likelihoods. Part of the issue of estimating huge trees is doing so in a reasonable timeframe. So if an analysis doesn't finish in 24 hours, I'll just can it as infeasible and note it below.

One last thing. From running replicate analyses on my server, I've noticed that I sometimes see a roughly 5% difference in runtime on replicate analyses. I have no idea quite why this occurs, but it does mean that we should probably ignore runtime differences that are lower than, say, about 10%.

I'll keep a note of the best settings in the table below. And I'll keep fully reproducible notes of the commandlines lower down.


## Results

### The best commandline so far

Here's the best commandline for inferring the ML tree. This only has to be done once, so we can use a lot of time and memory without worrying too much. 

```
fasttree -nosupport -nt -gamma global.fa
```



## Methods: model and topology search for fasttree

#### Iteration 1: Inititial analysis from 14/6/20

Add `-gamma` to make the likelihoods comparable between runs 

```
start=`date +%s`

fasttree -nosupport -nt -fastest -gamma global.fa > 1.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

#### Iteration 2

There are a couple of obvious ways I might be able to improve the tree. In this iteration I'll explore various combinations of using the `-fastest` flag, and using GTR or JC models. 

2.1: Don't use `-fastest`

```
start=`date +%s`

fasttree -nosupport -nt -gamma global.fa > 2.1.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

2.2: Use `-fastest` but with `-no2nd`, which makes the NJ phase a bit more rigorous

```
start=`date +%s`

fasttree -nosupport -nt -fastest -no2nd -gamma global.fa > 2.2.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

2.3: Use gtr instead of jc model

```
start=`date +%s`

fasttree -nosupport -gtr -nt -fastest -gamma global.fa > 2.3.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

Here are the GTR parameters that fasttree esimated:

```
ML-NNI round 1: LogLk = -398642.549 NNIs 14282 max delta 27.52 Time 20762.29es (max delta 27.520)      
GTR Frequencies: 0.2988 0.1836 0.1963 0.3213ep 12 of 12   
GTR rates(ac ag at cg ct gt) 0.2426 0.9987 0.1583 0.2036 2.9743 1.0000
Switched to using 20 rate categories (CAT approximation)20 of 20   
Rate categories were divided by 0.742 so that average rate = 1.0
CAT-based log-likelihoods may not be comparable across runs
  32550.38 seconds: ML NNI round 2 of 30, 1001 of 31287 splits, 198 changes (max delta 2.350)   
```



2.4 Use gtr with fixed parameters instead of jc model

We can do this with the matrix and base frequencies esimated by IQ-TREE. This may seem counterintuitive, but it's certainly better than assuming JC (equal weights, equal base frequencies). Also, on a 10K sequence alignment, the GTR parameters from IQ-TREE were:

```
Rate parameters:  A-C: 0.16442  A-G: 0.80624  A-T: 0.14639  C-G: 0.14328  C-T: 2.74138  G-T: 1.00000
Base frequencies:  A: 0.299  C: 0.183  G: 0.196  T: 0.322
Gamma shape alpha: 0.221
```

and on the much bigger and much later 40K sequence alignment from 14/6/20, the parameters were:

```
Rate parameters:  A-C: 0.17383  A-G: 0.75939  A-T: 0.16132  C-G: 0.14003  C-T: 2.80597  G-T: 1.00000
Base frequencies:  A: 0.299  C: 0.184  G: 0.196  T: 0.321
Gamma shape alpha: 0.372
```

So they really haven't changed a lot. Which means it's probably fine to fix them, and almost certainly better to fix them than to re-estimate them from scratch every update. Still, it's only worth fixing them if it actually saves time, so comparing the runtime of 2.6 to 2.3 will be useful.

```
start=`date +%s`

fasttree -nosupport -nt -fastest -gamma -gtrrates 0.17383 0.75939 0.16132 0.14003 2.80597 1.00000 -gtrfreq 0.299 0.184 0.196 0.321 global.fa > 2.4.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

###### Summary of iteration 2

It's pretty clear that turning off the -fastest option helped a lot (2.1) and that using a GTR model also helped a lot. The GTR model on its own did not fix the issue with the placement of lineage A.3 though, and it added about 20% to the runtime. Fixing the GTR parameters is a sensible way to address this, but oddly I noticed that for this run the memory usage ballooned out to about 3x the usage with the JT model. We'd want a reasonable improvement in the tree topology if we're spending that much extra RAM, so let's compare the topologies' likelihoods by estimating the lnL of the JC tree using a GTR model with the same parameters. Really, we want to know how much of the huge difference in likelihoods between analyses 2.1 and 2.3 is down to the tree, and how much to the model.

Let's do the analysis first:

```
fasttree -intree 2.1.tree -nome -mllen -nosupport -nt -gamma -gtrrates 0.2426 0.9987 0.1583 0.2036 2.9743 1.0000 -gtrfreq 0.2988 0.1836 0.1963 0.3213 global.fa > 2.5.tree
```

So the GTR tree (2.3) has a likelihood of -364798.195, while the JC tree (2.1) has a likelihood of -382048.463. That's a big difference of 17250. My suspicion is that a lot of that came from the GTR model parameters, and not much from the tree. Re-estimating the lnL of the JC tree with a GTR model gives a likelihood of -364581.620. This confirms something a bit surprising - the JC without the -fastest parameter tree is *better* than the GTR tree with it. So, as expected, the initial tree is making a lot of difference, and the ML NNIs at the end don't do so much. So, fine to keep the JC trees for bootstraps, but probably better to use the more rigorous options and the GTR model for the best tree (the branch lengths will be more accurate, for one). 

#### Iteration 3

Now we've tried some basic settings, it makes sense to see if we can further optimise things by making the initial topology better. In these analyses, `fasttree` typically spends about 20% of its runtime estimating the initial topology, and this is done with a fast version of Neighbour Joining. So if we're going to be able to save time, it's probably going to be at this step.

There are three parameters we can alter here. `-topm` (default 1.0) which affects the size of the list of top hits during NJ. `-close` (default 0.75) which is not explained much, but we are told that 'lower is more conservative'. `-refresh` (default 0.8, or 0.5 with `-fastest`), which affects which joined nodes are compared. It's really not clear to me how these three parameter will interact. Given that I already suspect that the tree is very sensitive to the value of the `-close` parameter, it seems silly to try to optimise the parameters one at a time. What I want to do ideally is find a combination that is accurate and fast. So instead I have to try combinations of parameter values to account for potential interactios. 

Let's start by seeing how low we can push the 'refresh' setting while maintaining reasonable likelihoods, and while making sure that lineage A.3 stays in the right place at the root of the tree.


3.1. refresh 0.8

```
start=`date +%s`

fasttree -refresh 0.8 -nosupport -nt -gamma global.fa > 3.1.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
nw_reroot 3.1.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 3.1.r.tree
```

3.2. refresh 0.7

```
start=`date +%s`

fasttree -refresh 0.7 -nosupport -nt -gamma global.fa > 3.2.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
nw_reroot 3.2.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 3.2.r.tree
```

3.3. refresh 0.6

```
start=`date +%s`

fasttree -refresh 0.6 -nosupport -nt -gamma global.fa > 3.3.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
nw_reroot 3.3.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 3.3.r.tree
```

3.4. refresh 0.5

```
start=`date +%s`

fasttree -refresh 0.5 -nosupport -nt -gamma global.fa > 3.4.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
nw_reroot 3.4.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 3.4.r.tree
```

3.5. refresh 0.4

```
start=`date +%s`

fasttree -refresh 0.4 -nosupport -nt -gamma global.fa > 3.5.tree

end=`date +%s`
runtime=$((end-start))
echo $runtime
nw_reroot 3.5.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 3.5.r.tree
```

###### Summary of iteration 3

Well, the timings are really odd here. I can only assume someone else snuck onto the server and started running lots of stuff. Oh well. Point is there are no really meaningful differences in timing or likelihood here. So it seems like the 'refresh' parameter is not doing much that matters for this tree. So I think setting it lower is proably better, since it's set lower with '-fastest' which I guess means it might save some memory. I'll track memory usage explicitly in future iteration to try and confirm or deny this.

#### Iteration 4

I've learned a new trick so I can track time and memory usage better, so I'll use that from here. In this iteration I'm just going to compare a few options that have occurred to me after reading and thinking more about how fasttree works.

For this iteration I'm going to switch back to using `-fastest` a bit, and see if I can get any improvements via some other methods, since switching off `-fastest` increases the execution time a bunch.

4.1 The best method so far (repeated to track memory usage)

```
/usr/bin/time -v fasttree -refresh 0.4 -nosupport -nt -gamma global.fa > 4.1.tree
nw_reroot 4.1.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 4.1.r.tree
```

4.2 Using a pre-set GTR
this gave better likelihood and was of course faster than letting fasttree estimate GTR parameters. I'm re-running it to double check the memory usage.

```
/usr/bin/time -v fasttree -refresh 0.4 -nosupport -nt -gamma -gamma -gtrrates 0.17383 0.75939 0.16132 0.14003 2.80597 1.00000 -gtrfreq 0.299 0.184 0.196 0.321 global.fa > 4.2.tree
nw_reroot 4.2.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 4.2.r.tree
```

4.3 Using a pre-set GTR, -fastest, and a couple more rounds of SPR moves

In principle SPR moves could fix issues with the misplacement of the A3 lineage. The default is 2. Let's increase this to 6, which is one round of SPR for every 10 NNI moves at the current dataset size. I'm also interested in the memory usage here, from memory it was the ML NNIs which cost a lot of memory here.

```
/usr/bin/time -v fasttree -fastest -refresh 0.4 -nosupport -nt -gamma -gamma -gtrrates 0.17383 0.75939 0.16132 0.14003 2.80597 1.00000 -gtrfreq 0.299 0.184 0.196 0.321 -spr 6 global.fa > 4.3.tree
nw_reroot 4.3.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 4.3.r.tree
```

4.4 Using a pre-set GTR, -fastest, and larger SPR moves of size 20

In principle SPR moves could fix issues with the misplacement of the A3 lineage. Instead of more of them, we could have the same number but make them longer (so subtrees can move further). The default is 10, so let's move to 20.

```
/usr/bin/time -v fasttree -fastest -refresh 0.4 -nosupport -nt -gamma -gamma -gtrrates 0.17383 0.75939 0.16132 0.14003 2.80597 1.00000 -gtrfreq 0.299 0.184 0.196 0.321 -sprlength 20 global.fa > 4.4.tree
nw_reroot 4.5.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 4.5.r.tree
```

4.5 Using a pre-set GTR, -fastest, and larger SPR moves of size 50

In principle SPR moves could fix issues with the misplacement of the A3 lineage. Instead of more of them, we could have the same number but make them longer (so subtrees can move further). The default is 10, so let's move to 20.

```
/usr/bin/time -v fasttree -fastest -refresh 0.4 -nosupport -nt -gamma -gamma -gtrrates 0.17383 0.75939 0.16132 0.14003 2.80597 1.00000 -gtrfreq 0.299 0.184 0.196 0.321 -sprlength 50 global.fa > 4.5.tree
nw_reroot 4.5.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 4.5.r.tree
```


###### Summary of iteration 4

Adding 4 more rounds of SPR moves increased the runtime a lot, and didn't help the likelihood at all (it made it a bit worse). But increasing the depth of the 2 SPR moves from 10 to 20 did help - improved the runtime and didn't affect the lnL. So, next up I'll look in more detail at that parameter, and switch back to the JC model because the GTR model uses so much RAM (45GB per run at peak RAM usage).



#### Iteration 5

Now I'm just going to try a few combinations of SPR depth (10, 15, 20, 25, 30, 35, 40) and the `-refresh` parameter (0.4, 0.2, 0.1). That's a lot of analyses (21), but it should be informative. Here are all the commandlines, which for simplicity I will just bung in a text file called `iteration5.txt` and run with:

```
parallel -j 21 -k --bar :::: iteration5.txt
```

this just runs the jobs in blocks of 21, keeping them in the same order as below, and printing a progress bar.

```
/usr/bin/time -v fasttree -sprlength 10 -refresh 0.4 -nosupport -nt -gamma -log 5.1.txt global.fa > 5.1.tree
/usr/bin/time -v fasttree -sprlength 15 -refresh 0.4 -nosupport -nt -gamma -log 5.2.txt global.fa > 5.2.tree
/usr/bin/time -v fasttree -sprlength 20 -refresh 0.4 -nosupport -nt -gamma -log 5.3.txt global.fa > 5.3.tree
/usr/bin/time -v fasttree -sprlength 25 -refresh 0.4 -nosupport -nt -gamma -log 5.4.txt global.fa > 5.4.tree
/usr/bin/time -v fasttree -sprlength 30 -refresh 0.4 -nosupport -nt -gamma -log 5.5.txt global.fa > 5.5.tree
/usr/bin/time -v fasttree -sprlength 35 -refresh 0.4 -nosupport -nt -gamma -log 5.6.txt global.fa > 5.6.tree
/usr/bin/time -v fasttree -sprlength 40 -refresh 0.4 -nosupport -nt -gamma -log 5.7.txt global.fa > 5.7.tree
/usr/bin/time -v fasttree -sprlength 10 -refresh 0.2 -nosupport -nt -gamma -log 5.8.txt global.fa > 5.8.tree
/usr/bin/time -v fasttree -sprlength 15 -refresh 0.2 -nosupport -nt -gamma -log 5.9.txt global.fa > 5.9.tree
/usr/bin/time -v fasttree -sprlength 20 -refresh 0.2 -nosupport -nt -gamma -log 5.10.txt global.fa > 5.10.tree
/usr/bin/time -v fasttree -sprlength 25 -refresh 0.2 -nosupport -nt -gamma -log 5.11.txt global.fa > 5.11.tree
/usr/bin/time -v fasttree -sprlength 30 -refresh 0.2 -nosupport -nt -gamma -log 5.12.txt global.fa > 5.12.tree
/usr/bin/time -v fasttree -sprlength 35 -refresh 0.2 -nosupport -nt -gamma -log 5.13.txt global.fa > 5.13.tree
/usr/bin/time -v fasttree -sprlength 40 -refresh 0.2 -nosupport -nt -gamma -log 5.14.txt global.fa > 5.14.tree
/usr/bin/time -v fasttree -sprlength 10 -refresh 0.1 -nosupport -nt -gamma -log 5.15.txt global.fa > 5.15.tree
/usr/bin/time -v fasttree -sprlength 15 -refresh 0.1 -nosupport -nt -gamma -log 5.16.txt global.fa > 5.16.tree
/usr/bin/time -v fasttree -sprlength 20 -refresh 0.1 -nosupport -nt -gamma -log 5.17.txt global.fa > 5.17.tree
/usr/bin/time -v fasttree -sprlength 25 -refresh 0.1 -nosupport -nt -gamma -log 5.18.txt global.fa > 5.18.tree
/usr/bin/time -v fasttree -sprlength 30 -refresh 0.1 -nosupport -nt -gamma -log 5.19.txt global.fa > 5.19.tree
/usr/bin/time -v fasttree -sprlength 35 -refresh 0.1 -nosupport -nt -gamma -log 5.20.txt global.fa > 5.20.tree
/usr/bin/time -v fasttree -sprlength 40 -refresh 0.1 -nosupport -nt -gamma -log 5.21.txt global.fa > 5.21.tree
```

Then I'll re-root the trees similarly:

```
parallel -j 21 -k --bar :::: iteration5root.txt
```

with this text file:

```
nw_reroot 5.1.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.1.r.tree
nw_reroot 5.2.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.2.r.tree
nw_reroot 5.3.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.3.r.tree
nw_reroot 5.4.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.4.r.tree
nw_reroot 5.5.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.5.r.tree
nw_reroot 5.6.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.6.r.tree
nw_reroot 5.7.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.7.r.tree
nw_reroot 5.8.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.8.r.tree
nw_reroot 5.9.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.9.r.tree
nw_reroot 5.10.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.10.r.tree
nw_reroot 5.11.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.11.r.tree
nw_reroot 5.12.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.12.r.tree
nw_reroot 5.13.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.13.r.tree
nw_reroot 5.14.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.14.r.tree
nw_reroot 5.15.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.15.r.tree
nw_reroot 5.16.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.16.r.tree
nw_reroot 5.17.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.17.r.tree
nw_reroot 5.18.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.18.r.tree
nw_reroot 5.19.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.19.r.tree
nw_reroot 5.20.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.20.r.tree
nw_reroot 5.21.tree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > 5.21.r.tree
```




#### Iteration 6

The results on the placement of clade A3 are still baffling. More and/or bigger SPR moves in iteration 4 puts it back in the wrong spot (confirmed probably by the trees having lower likelihoods), and that's a fairly big worry. But my suspicion is that this is driven by the GTR model placing high weight on crazy hypermutable sites. Switching back to the JC model addresses this to some extent. Iteration 5 suggests that there's not a super strong effect of increasing the SPR length, but this doesn't seem to have a huge effect on timing either. So it seems smarter to keep the SPR length high. The best setting was to have it at 20 - this gave good results regardless of the refresh setting.

So, now we try to improve the NJ tree without blowing the RAM usage out of the water. Since a large proportion of the execution time is spent on the NJ tree, and since all other steps have to work to fix issues introduced in this tree, making it better is probably a good investment in accuracy and potentially in time.

Now I'll try: 

* `refresh` of: 0.8, 1.0, 1.2
* `topm` of: 1.0, 1.2, 1.5
* `close` of 0.75 and 0.5

These are the three parameters I have access to that ajust the effor that goes into the NJ tree. I don't know how they interact, how they affect execution time, or how they affect memory. Neither do I know whether making the NJ tree a better NJ tree (which is what these can do) makes it a better ML tree. These analyses will help me figure some of that out. 

```
/usr/bin/time -o 6.1.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.0 -close 0.75 -nosupport -nt -gamma -log 6.1.txt global.fa > 6.1.tree
/usr/bin/time -o 6.2.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.0 -close 0.5 -nosupport -nt -gamma -log 6.2.txt global.fa > 6.2.tree
/usr/bin/time -o 6.3.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.2 -close 0.75 -nosupport -nt -gamma -log 6.3.txt global.fa > 6.3.tree
/usr/bin/time -o 6.4.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.2 -close 0.5 -nosupport -nt -gamma -log 6.4.txt global.fa > 6.4.tree
/usr/bin/time -o 6.5.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.75 -nosupport -nt -gamma -log 6.5.txt global.fa > 6.5.tree
/usr/bin/time -o 6.6.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.5 -nosupport -nt -gamma -log 6.6.txt global.fa > 6.6.tree
/usr/bin/time -o 6.7.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.0 -close 0.75 -nosupport -nt -gamma -log 6.7.txt global.fa > 6.7.tree
/usr/bin/time -o 6.8.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.0 -close 0.5 -nosupport -nt -gamma -log 6.8.txt global.fa > 6.8.tree
/usr/bin/time -o 6.9.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.2 -close 0.75 -nosupport -nt -gamma -log 6.9.txt global.fa > 6.9.tree
/usr/bin/time -o 6.10.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.2 -close 0.5 -nosupport -nt -gamma -log 6.10.txt global.fa > 6.10.tree
/usr/bin/time -o 6.11.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.5 -close 0.75 -nosupport -nt -gamma -log 6.11.txt global.fa > 6.11.tree
/usr/bin/time -o 6.12.mem.txt -v fasttree -sprlength 20 -refresh 0.9 -topm 1.5 -close 0.5 -nosupport -nt -gamma -log 6.12.txt global.fa > 6.12.tree
/usr/bin/time -o 6.13.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.0 -close 0.75 -nosupport -nt -gamma -log 6.13.txt global.fa > 6.13.tree
/usr/bin/time -o 6.14.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.0 -close 0.5 -nosupport -nt -gamma -log 6.14.txt global.fa > 6.14.tree
/usr/bin/time -o 6.15.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.2 -close 0.75 -nosupport -nt -gamma -log 6.15.txt global.fa > 6.15.tree
/usr/bin/time -o 6.16.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.2 -close 0.5 -nosupport -nt -gamma -log 6.16.txt global.fa > 6.16.tree
/usr/bin/time -o 6.17.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.5 -close 0.75 -nosupport -nt -gamma -log 6.17.txt global.fa > 6.17.tree
/usr/bin/time -o 6.18.mem.txt -v fasttree -sprlength 20 -refresh 0.999 -topm 1.5 -close 0.5 -nosupport -nt -gamma -log 6.18.txt global.fa > 6.18.tree
```


The upshot of this is that 6.5 gives the best results, with refresh of 0.9, topm of 1.5 and close of 0.75. 


#### Iteration 7

Let's try a few more values of topm, since we didn't reach a maximum value of that parameter, and increasing it didn't increase the execution time appreciably.

* topm: 1.6, 1.8, 2.0, 2.25, 2.5

```
/usr/bin/time -o 7.1.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.6 -close 0.75 -nosupport -nt -gamma global.fa > 7.1.tree
/usr/bin/time -o 7.2.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 1.8 -close 0.75 -nosupport -nt -gamma global.fa > 7.2.tree
/usr/bin/time -o 7.3.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.0 -close 0.75 -nosupport -nt -gamma global.fa > 7.3.tree
/usr/bin/time -o 7.4.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -nosupport -nt -gamma global.fa > 7.4.tree
/usr/bin/time -o 7.5.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.5 -close 0.75 -nosupport -nt -gamma global.fa > 7.5.tree
```

The best settings are for 7.4, which is quick and gets a better likelihood.

This means that the NJ tree is now fairly well optimised. So let's move on to the ME tree after trying a few more options that I'd like to see on the NJ tree.

#### Iteration 8

A few more random things to try on the NJ tree. Mostly because I'm curious.

```
/usr/bin/time -o 8.1.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -nt -gamma global.fa > 8.1.tree
/usr/bin/time -o 8.2.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -bionj -nosupport -nt -gamma global.fa > 8.2.tree
/usr/bin/time -o 8.3.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -2nd -nosupport -nt -gamma global.fa > 8.3.tree
/usr/bin/time -o 8.4.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -bionj -2nd -nosupport -nt -gamma global.fa > 8.4.tree
/usr/bin/time -o 8.5.mem.txt -v fasttree -sprlength 20 -refresh 0.8 -topm 2.25 -close 0.75 -fastest -bionj -nosupport -nt -gamma global.fa > 8.5.tree
```

Nothing there helped!


#### Iteration 9

Big switch here. For various reasons I decided to switch to only using alignments from GISAID that are 'complete' and 'high coverage', according to GISAID's guidlines. So I'm switching the alignment to one of those from 24/6/20. That means likelihoods from iteration 9 are not comparable to those from previous iterations.

Let's start by estimating the tree using our best efforts from previous iterations.

```
/usr/bin/time -o 9.0.mem.txt -v fasttree -nt -gamma -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.75 global.fa > 9.0.tree
```

The likelihood of this tree is -330618.967

Now I'm intersted to see if I can improve on that tree with some additional SPR moves, but this will require a little bit of hacking around with fasttree.

Here's the idea. How well can we do with three rounds of fasttree, starting with two rounds of big SPR moves on the ME tree, and finishing with a standard analysis. 

```
/usr/bin/time -o 9.1.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 100 -nni 0 -spr 10 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 9.1.1.tree
/usr/bin/time -o 9.1.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 50  -nni 0 -spr 10 -refresh 0.8 -topm 1.5 -close 0.75 -noml -intree 9.1.1.tree global.fa > 9.1.2.tree
/usr/bin/time -o 9.1.3.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.75 -intree 9.1.2.tree global.fa > 9.1.3.tree
```


Would it be better to start with the ML tree 9.0 instead?

```
/usr/bin/time -o 9.2.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 100 -nni 0 -spr 10 -refresh 0.8 -topm 1.5 -close 0.75 -noml -intree 9.0.tree global.fa > 9.2.1.tree
/usr/bin/time -o 9.2.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 50  -nni 0 -spr 10 -refresh 0.8 -topm 1.5 -close 0.75 -noml -intree 9.2.1.tree global.fa > 9.2.2.tree
/usr/bin/time -o 9.2.3.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.75 -intree 9.2.2.tree global.fa > 9.2.3.tree
```

The first option (starting from scratch) gives us a likelihood of -330458.212, which is a decent improvement. The second option gives us -330481.132, which is similar but not quite as good.

The first option took 36259.18 seconds for the first step (long because this includes the NJ tree inference, which is slow), 13896.96 for the second step, and 13511.66 for the third step. This is fine for the ML tree, the total time is blown out a bit, but these long-range SPR moves should help avoid potentially large issues in the resulting tree. In short, additional SPR rounds do not take much time, and are certainly worth it for reducing errors.

So the upshot is - extra SPR moves are a good thing. Now let's see just how much difference it made at each step, to guide our intuition about how many SPR moves we actually need.

```
fasttree -nt -gamma -nosupport -mllen -nome -intree 9.0.tree global.fa > ignore1.tree
fasttree -nt -gamma -nosupport -mllen -nome -intree 9.1.1.tree global.fa > ignore2.tree
fasttree -nt -gamma -nosupport -mllen -nome -intree 9.1.2.tree global.fa > ignore3.tree
fasttree -nt -gamma -nosupport -mllen -nome -intree 9.1.3.tree global.fa > ignore4.tree
```


#### Iteration 10

So now we know that a few extra rounds of SPR moves are a good idea, we need to figure out just how many SPR moves are worth doing, and at what range.

An obvious question is whether it's just always better to have lots of very long range SPRs, so let's try that. Indeed, the range of the SPR moves didn't seem to change the execution time much, so it's feasible that increasing the range won't slow us down. Let's stick with 20 SPR moves before the final ML infernce, but adjust the range to a few different levels.

And from a gut feeling, it also seems likely to improve things if we intersperse a few NNI moves between the SPR moves, so let's add two rounds of NNI moves for every round of SPR moves.

Once we find the best of these options, it might be worth trying to save some execution time by reducing the number of rounds of SPR and NNI moves. 

```
/usr/bin/time -o 10.1.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 500 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.1.1.tree
/usr/bin/time -o 10.1.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 500 -intree 10.1.1.tree global.fa > 10.1.2.tree
```

```
/usr/bin/time -o 10.2.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 200 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.2.1.tree
/usr/bin/time -o 10.2.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 200 -intree 10.2.1.tree global.fa > 10.2.2.tree
```

```
/usr/bin/time -o 10.3.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 100 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.3.1.tree
/usr/bin/time -o 10.3.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 100 -intree 10.3.1.tree global.fa > 10.3.2.tree
```

```
/usr/bin/time -o 10.4.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 50 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.4.1.tree
/usr/bin/time -o 10.4.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 50 -intree 10.4.1.tree global.fa > 10.4.2.tree
```

```
/usr/bin/time -o 10.5.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 20 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.5.1.tree
/usr/bin/time -o 10.5.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 20 -intree 10.5.1.tree global.fa > 10.5.2.tree
```

```
/usr/bin/time -o 10.6.1.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 10 -nni 40 -spr 20 -refresh 0.8 -topm 1.5 -close 0.75 -noml global.fa > 10.6.1.tree
/usr/bin/time -o 10.6.2.mem.txt -v fasttree -nt -gamma -nosupport -sprlength 10 -intree 10.6.1.tree global.fa > 10.6.2.tree
```
