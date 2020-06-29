# How fast and accurate is VeryFastTree compared to fasttree?

### TL;DR

VeryFastTree *is* faster than single-threaded fasttree on GISAID alignments, by aout 40% if you use 5 threads. But the trees are worse and it uses 3 times the memory. 

Also, you can have fun guessing what the next attempt to speed this program up will be called. Our lab Slack has been proposing things like `UltraFastTree`, `TolerablySlowTree` (this is not meant to be a dig at the software, more at naming conventions of which we are certainly guilty having between us published things like the UltraFast Bootstrap and Fast Tiger...), etc. 

### What's this?

With impeccable timing, a group just released a new version of fasttree called VeryFastTree, which purports to be identical to fasttree but much faster in a few key ways, by using vectorisation and better multithreading.

Let's try it out. I'll run it with all the obvious different options, time it (roughly - this is on a server that is also used by others so the timings are not going to be perfectly comparable), look at the trees we get and the memory used.

Happily, it compiled first time, with no issues at all. That's a good start. Oddly though when I ran it first time I tried to specify the `close` parameter to be `0.75` I got the following nonsensical error:

`-close: -close argument must be between 0 and 1`

That's a bit of a worry - if the software has straightforward bugs in input parameters I find it a little harder to trust that it's going to work as advertised. Regardless - onwards and upwards! Let's see if we can get it to run at all. Turns out if you don't use any of the options for controlling the NJ aspect of treebuilding, it runs OK.

So, let's compare it to fasttree. I'll jump straight to the results. Methods are further down.

### Results

| program      | threads | ext  | fastexp | runtime  | max.memory | lnL         | bad.splits |
|--------------|---------|------|---------|----------|------------|-------------|------------|
| fasttree     | 1       | NA   | NA      | 46986.38 | 13034040   | -330691.544 | 0          |
| veryfasttree | 1       | AVX2 | 0       | 45651.11 | 36281580   | -330764.426 | 111        |
| veryfasttree | 2       | AVX2 | 0       | 33248.28 | 36300084   | -330764.426 | 111        |
| veryfasttree | 3       | AVX2 | 0       | 29159.18 | 36518228   | -330755.747 | 108        |
| veryfasttree | 4       | AVX2 | 0       | 29498.8  | 36554492   | -330726.65  | 105        |
| veryfasttree | 5       | AVX2 | 0       | 27357.99 | 36730576   | -330882.389 | 123        |
| veryfasttree | 1       | SSE3 | 0       | 46817.54 | 36290956   | -330758.843 | 112        |
| veryfasttree | 1       | SSE3 | 2       | 47008.52 | 36281572   | -330764.426 | 111        |
| veryfasttree | 1       | AVX  | 0       | 43583.64 | 36296364   | -330758.843 | 112        |
| veryfasttree | 1       | AVX  | 2       | 46528.24 | 36296412   | -330758.843 | 112        |
| veryfasttree | 1       | AVX2 | 0       | 47030.92 | 36296484   | -330758.843 | 112        |
| veryfasttree | 1       | AVX2 | 2       | 46741.94 | 36296424   | -330758.843 | 112        |


### Conclusions

VeryFastTree is faster, by a little bit. Note that a really good comparison would compare VeryFastTree to the multi-threaded fasttree. I didn't do that. I compared it only to the single-threaded fasttree. Nevertheless, if you give VeryFastTree 5 threads, it's shaves about 40% off the running time compared to single-threaded fasttree. That's OK but not stunning. It's a roughly 1.7x speedup. Perhaps the runtime benefits of VeryFastTree are only really apparent on HUGE alignments - the paper abstract mentions a 7.8x speedup on a 330k species alignment (this one is just 35k). 

Notably (perhaps) forcing VeryFastTree to use different vectorisations (the `ext` parameter) made no difference to the runtime (remember the timings are rough, so probably only accurate within a couple of thousand seconds), and neither did using the supposedly faster matrix exponentiation trick (`fastexp` of 2). 

VeryFastTree always uses 3x more memory than fasttree. That's a problem for me since memory is limiting when I run 100 true bootstraps.

Most worryingly for me, VeryFastTree always finds worse trees than fasttree. If I've understood correctly, the likelihoods should be directly comparable between the two programs. If that's true, then VeryFastTree was always worse than fasttree. Most worryingly, the very worst tree (no pun intended) came with the fastest setting using 5 threads. Since I really care about the quality of the tree, this is a genuine concern for me. VeryFastTree also returns trees with a few bad splits (I'm guessing these are splits where a NNI move that's *not* in the tree is preferred to one that is). I'm not sure why this is, but it seems to just confirm that the VeryFastTree trees are not as well optimised as the fasttree trees.


### Methods

These comparisons are all done on all of the complete and high coverage genomes from GISAID on the 24th of June 2020. I can't share the alignment publicly, but if you're a GISAID member you can email me and I'll send it to you. Otherwise you can build it by following the steps outlined in this repo. 

First, we'll fun fasttree as normal. 

```
/usr/bin/time -o fasttree.mem.txt -v fasttree -nt -gamma global.fa > fasttree.tree
```

Now we'll try VeryFastTree in a few different ways. Very nicely, the authors made it so that the same commandline should just work. So let's start there. Note that I don't need to specify the number of threads, but I will just to be sure I know what's happening. However, I do need to explicitly tell VeryFastTree to use double precision (crucial for SARS-CoV-2 alignment because they have so little variation and such short branches)

```
/usr/bin/time -o veryfasttree1.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision global.fa > veryfasttree1.tree 
```

But there are a few other options we can also try, listed here: https://github.com/citiususc/veryfasttree#running-veryfasttree

Let's see how it goes as we increase the number of threads from 1-5 (realistically, this won't help my cause much, because I need to estimate 100 bootstrap trees on a server with 56 threads, but it's still nice to know how it goes).

```
/usr/bin/time -o veryfasttree2.mem.txt -v VeryFastTree -nt -gamma -threads 2 -double-precision global.fa > veryfasttree2.tree 

/usr/bin/time -o veryfasttree3.mem.txt -v VeryFastTree -nt -gamma -threads 3 -double-precision global.fa > veryfasttree3.tree 

/usr/bin/time -o veryfasttree4.mem.txt -v VeryFastTree -nt -gamma -threads 4 -double-precision global.fa > veryfasttree4.tree 

/usr/bin/time -o veryfasttree5.mem.txt -v VeryFastTree -nt -gamma -threads 5 -double-precision global.fa > veryfasttree5.tree 
```

And let's try a few different architectures (my server has AVX2 chips), combined with two different options for speeding up the exponentiation (only options 0 and 2 are feasible with double precision, according to the documentation). We'll do this all on one thread for simplicity.

```
/usr/bin/time -o veryfasttree6.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext SSE3 -fastexp 0 global.fa > veryfasttree6.tree 

/usr/bin/time -o veryfasttree7.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext SSE3 -fastexp 2 global.fa > veryfasttree7.tree 

/usr/bin/time -o veryfasttree8.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext AVX -fastexp 0 global.fa > veryfasttree8.tree 

/usr/bin/time -o veryfasttree9.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext AVX -fastexp 2 global.fa > veryfasttree9.tree 

/usr/bin/time -o veryfasttree10.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext AVX2 -fastexp 0 global.fa > veryfasttree10.tree 

/usr/bin/time -o veryfasttree11.mem.txt -v VeryFastTree -nt -gamma -threads 1 -double-precision -ext AVX2 -fastexp 2 global.fa > veryfasttree11.tree 
```