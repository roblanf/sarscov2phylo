# IQ-TREE of fasttree?

In the first benchmark of methods for SARS-CoV-2 phylogenetic analyses, fasttree was the clear winner, and I've been using it succesfully for months. But IQ-TREE is certianly a more flexible piece of software, and what's more we develop it. As luck would have it, we have funding from the Chan-Zuckerberg Initiative, and we used it to hire an *excellent* programmer, James Barbetti. The first thing James has been looking at is removing some of the bottlenecks that IQ-TREE had when analysing huge datasets. 

The result so far is that we have a version of IQ-TREE that does lightning fast neighbour joining trees. This is done using all sorts of methods, some of which come straight from the excellent work that went into developing the rapid-NJ algorithm used in `rapidnj`. Regardless, now's a good time to revisit whether we can get better trees from fasttree or IQ-TREE. Of note, I've spent month optimising the search in fasttree (see other .md files on this repo), to the point where I don't think I can make much in the way of additional gains. 

So I'm not expecting IQ-TREE to beat fasttree first time. It will likely take a lot of tweaks. It may not be possible at all, because fasttree still does things IQ-TREE can't, like SPR moves, which at least in fasttree seem to be fairly crucial to getting good and accurate trees from SARS-CoV-2. But I still want to check.

The other reason for this is that IQ-TREE allows us to apply a much much wider range of models of sequence evolution, which I know from previous analyses fit the data way better than the JC+G or GTR+G models availble in fasttree. In addition, we are developing a range of new models in IQ-TREE, specifically targetted at dealing with the various issues we see in SARS-CoV-2 sequencing data (see posts on virological that highlight these).

## The plan

I'll start by figuring out a good model of sequence evolution in IQ-TREE. Based on previous analyses, the best I've found is `JC+I+R3`. But this was on an old alignment with different filtering, so it needs to be revisited. And I didn't even check the JC component of this. 

A few open questions:

1. What's the best rate matrix framework to use (e.g. JC, HKY, GTR, UNREST)?
2. What's the best distribution of rates across sites (e.g. +I+R2, +I+R4 etc)?
3. Do either 1 or 2 make much difference to the tree topology? (this is important - it's much much harder to optimise the complex models, so it may be better to spend a long time optimising the tree on a simple model like JC first, then do a few last rounds of optimisation with a complex model)
4. What tree search settings are best in IQ-TREE?
4. Does the best tree from IQ-TREE beat the best tree from fasttree (when both are analysed in IQ-TREE under whatever the best model seems to be).


These questions suggest a simple approach:

1. Find the best rate distribution
2. Find the best rate matrix
3. Optimise in IQ-TREE with this model, and with plain JC model
3. Get the lnL of the best fasttree tree under this model
5. Compare fasttree to the IQ-TREE trees
6. Work on optimising IQ-TREE

Doing things in this order makes sense. The finder details of the tree topology tend not to make much difference for the model selection steps. And getting the rate distribution right tends to be more important than getting the rate matrix right. I'm not going to worry about running time here, because ultimately I can fix all the model parameters in IQ-TREE, I know they don't change much as the size of the dataset increases, so in practice my analyses should be much faster than the analyses here. In other words, at this stage I care about lnL values much more than I care about run time. 

Of note, IQ-TREE uses a lot more memory than fasttree. So I'm very limited in how many options I can explore here. Typically I can just look at 4 options at a time. With that plain laid out, here goes. 

I'll do all of the analyses on the latest global alignment from this analysis: https://github.com/roblanf/sarscov2phylo/releases/tag/11-7-20

Also of note, I am using here the bleeding edge IQ-TREE release with all sorts of udocumented options, which is version 2.0.8. For all analyses I'll set the minimum branch length to a very small number.

## 1. Find the best rate distribution

I know from previous analyses that the +G and +I+G models are very unlikely to be the best. So I'll just compare four models here that I think should capture a good range.

```
/usr/bin/time -o 1.1.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.1

/usr/bin/time -o 1.2.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I+R3 -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.2

/usr/bin/time -o 1.3.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I+R4 -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.3

/usr/bin/time -o 1.4.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I+R5 -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.4
```

##### Results

Don't pay too much heed to execution time. We're working on speeding up bottlenecks. 

JC 
lnL -407738.7183
BIC: 1464697.832
time: 6h:53m:34s

JC+I+R3
lnL: -388587.635
BIC: 1426447.16
time: 53h:37m:6s

JC+I+R4
lnL: -388425.334
BIC: 1426143.156
time: 105h:52m:55s

JC+I+R5
lnL: ? it was converging on something a little less than -388319.099658 though
BIC: ?
time: >200h (!)

Although the final run didn't complete (I had to kill it to re-estimate the global tree, which was a week out of date...), it did reach a lnL of -388319.099658, which equates to a BIC that is also substantially better than the JC+I+R4 model. 

It's worth comparing the model parameters estimated by IQ-TREE (the +R5 model parameters never completed their final optimisation, but even before that would have had the best BIC score) just to see what's going on.

```
Proportion of invariable sites: 0.389 ;Site proportion and rates:  (0.529,0.808) (0.074,5.739) (0.009,16.933)
Proportion of invariable sites: 0.335; Site proportion and rates:  (0.558,0.659) (0.083,4.126) (0.023,10.607) (0.002,27.445)
Proportion of invariable sites: 0.305 Site proportion and rates:  (0.366,0.494) (0.220,0.840) (0.084,4.024) (0.024,10.400) (0.002,27.193)
```

The rate categories show that there are roughly 30% of sites modelled as invariant. Then there's about 30% of sites rate around 0.5, 20% at 0.8, 8% at 4, 2% at 10, and 0.2% at 27%.

That last category represents around 60 sites that seem to have exceptionally fast rates relative to all the other sites. Modelling those sites properly (i.e. with a free rate model) could be important for correctly placing sequences in the tree.

Given that the I+R5 model gave the best BIC score, it seems prudent to use this model, but to fix the model parameters to those estimated here. They won't be the ML parameters (because they never quite finished the final optimisation, and fixing them means they won't be updated when the tree changes), however, it's a safe bet from the BIC scores that they are better than the other models I tried, and a *lot* better than the JC model with no rates. This could end up being important for optimising the tree - as I currently use fasttree, which doesn't account for rate variation when optimising the tree (it adds it at the end when optimising branch lengths).

Final note: the memory usage for the +I+R5 model is HUGE. 175GB at the maximum (which was most of the time)!

Postscript: turns out IQ-TREE has a potentially much faster algorithm for the model optimisation implemented. So let's try and compare on the +R3 and +R5 models (the former so I can compare likelihoods, the latter to see if it's quick enough to do this). The standard algorithm tends to perform better on most datasets, but it's slower. Usually that's not an issue. The standard algorithm updates the rates one at a time, jiggling that rate parameter around to improve it while holding the others constant. The alternative BFGS algorithm tries out changes to all rate parameters at once. Let's see what happens...


```
/usr/bin/time -o 1.5.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I+R3 -fast -nt 4 -optalg 1-BFGS -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.5

/usr/bin/time -o 1.6.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I+R5 -fast -nt 4 -optalg 1-BFGS -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.6

```

OK well that was a surprise. My expectation was that the BFGS algorithm would be quicker, but get substantially worse likelihoods. The truth was the opposite. I had to kill the runs because they took a long time (>> 1 day, so I don't know the exact speed comparison but the BFGS is certainly not more than twice as fast, and from looking at the stage it got to, I'd estimate it's about the same speed in fact). Speed aside, one can see the likelihoods are a lot better, even when it was just the model parameters optimised on the initial JC NJ-R tree.  

JC+I+R3:  -375974.923841
JC+I+R5:  -375490.305847

The difference between R3 and R5 is large, plenty to justify the two extra parameters I'd argue. The main take home is to compare these likelihoods to those above. They are A LOT better. Like more than 10K likelihood units. For the SAME model. That's huge. 

The model parameters:

Proportion of invariable sites: 0.616; Site proportion and rates:  (0.426,0.269) (0.415,0.625) (0.159,3.936)
Proportion of invariable sites: 0.616; Site proportion and rates:  (0.210,0.334) (0.219,0.342) (0.284,0.344) (0.202,1.456) (0.085,5.456)


These are very different to using the standard algorithm, and a lot better too (the likelihood tells us that). 

Then I realised I forgot one really important model. The plain-old +I model. Not sure quite how I forgot this, but let's see how it goes. We don't need BFGS for this. Won't help. 


```
/usr/bin/time -o 1.8.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m JC+I -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 1.8
```

JC+I: -396992.024799

The likelihood is a lot worse. But it did run quickly (about 5 hours). This may be worth considering in future.

So, let's continue with the +I+R5 model.


## 2. Find the best rate matrix

We can fix the rate parameters to the ones estimated in step 1 using:

`+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}`


I'll re-run the JC model too, to get an idea of how much time and memory we save by running this. 

```
/usr/bin/time -o 2.1.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "JC+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 2.1

/usr/bin/time -o 2.2.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "HKY+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 2.2

/usr/bin/time -o 2.3.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "GTR+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 2.3

/usr/bin/time -o 2.4.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "12.12+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -fast -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 2.4

```


lnLs (after first round of model optimisation on NJ-R ML distance tree)

JC+I+R5: -392823.153805
HKY+I+R5: -380249.398257
GTR+I+R5: -373470.588371
12.12+I+R%: didn't optimise, various numerical issues

So the GTR model is by far and away the best model. The parameters are: 

Rate parameters:  A-C: 0.14012  A-G: 0.66091  A-T: 0.12269  C-G: 0.09998  C-T: 2.61404  G-T: 1.00000
Base frequencies:  A: 0.299  C: 0.184  G: 0.196  T: 0.321
Proportion of invariable sites: 0.616
Site proportion and rates:  (0.077,0.280) (0.077,0.938) (0.077,1.823) (0.077,3.185) (0.077,6.795)


## 3. Optimise in IQ-TREE

Here I want to optimise the tree with the two different models. 

the -fast flag does this:

    #     -ninit 2
    #     -n 2
    #     -me 0.05

But I want to allow for more rounds of nni than this, so I'll allow for 5 rounds and hard-code the other options. 

I'll also set -me to 1.0, because I'm not really worried about optimising branch lengths at this point (and hte rest of the model parameters are fixed).

```
/usr/bin/time -o 3.1.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "JC+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -ninit 2 -n 5 -me 1.0 -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 3.1

/usr/bin/time -o 3.2.mem.txt -v ./iqtree-2.0.8-Linux/bin/iqtree2 -s global.fa -m "GTR{0.14012,0.66091,0.12269,0.09998,2.61404}+F{0.299,0.184,0.196,0.321}+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" -ninit 2 -n 5 -me 1.0 -nt 4 -experimental --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -t NJ-R --no-opt-gamma-inv -blmin 0.0000000001 -pre 3.2
```



