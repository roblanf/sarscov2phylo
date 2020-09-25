# Benchmarking fasttree OMP and veryfasttree (again)

Question - how much faster is fasttree OMP (compiled on my machine) than the non-OMP version I get from conda, and how many threads is optimal. And how does it compare to veryfasttree?

I looked at veryfasttree before (see veryfasttree.md on this repo), but at that point my trees were smaller, and I was limited by memory and CPUs, since I was running full bootstraps. Now my trees are much bigger, and I'm no longer doing full bootstraps, so it seems worth revisiting that analysis. 

# Methods

This one's simple. First I get fasttree OMP and veryfasttree and compile them:

```
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DOPENMP -fopenmp -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm

git clone https://github.com/citiususc/veryfasttree
cmake .
make
make install

```

Now I re-run the analyses from 14-09-20 and time them, setting different numbers of threads like this:

```
export OMP_NUM_THREADS=1
/usr/bin/time -o 1threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 1threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 1thread.tree
```

I'll just do this for 1, 2, 5, 10, 20, 50, and 100 threads like this:

```
export OMP_NUM_THREADS=2
/usr/bin/time -o 2threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 2threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 2thread.tree

export OMP_NUM_THREADS=5
/usr/bin/time -o 5threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 5threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 5thread.tree

export OMP_NUM_THREADS=10
/usr/bin/time -o 10threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 10threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 10thread.tree

export OMP_NUM_THREADS=20
/usr/bin/time -o 20threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 20threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 20thread.tree

export OMP_NUM_THREADS=50
/usr/bin/time -o 50threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 50threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 50thread.tree

export OMP_NUM_THREADS=100
/usr/bin/time -o 100threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 100threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 100thread.tree
```

I noticed during these runs that 5 threads seemed to be going quickest, so I thought I'd also fire up everything from 2-10, because if you're going to optimise, you may as well do it right.


```
export OMP_NUM_THREADS=3
/usr/bin/time -o 3threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 3threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 3thread.tree

export OMP_NUM_THREADS=4
/usr/bin/time -o 4threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 4threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 4thread.tree

export OMP_NUM_THREADS=6
/usr/bin/time -o 6threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 6threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 6thread.tree

export OMP_NUM_THREADS=7
/usr/bin/time -o 7threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 7threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 7thread.tree

export OMP_NUM_THREADS=8
/usr/bin/time -o 8threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 8threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 8thread.tree

export OMP_NUM_THREADS=9
/usr/bin/time -o 9threads.txt -v ../FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log 9threads.log -intree iqtree_seqsadded_mp.treefile global.fa > 9thread.tree
```

Now I'll run the same kinds of analysis for veryfasttree. But before I do that I'll run some analyses to find out which of veryfasttree's tricks work best on my new server:

```
/usr/bin/time -o vf1.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf1.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext SSE3 -fastexp 0 global.fa > vf1.tree 
/usr/bin/time -o vf2.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf2.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext SSE3 -fastexp 2 global.fa > vf2.tree 
/usr/bin/time -o vf3.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf3.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext AVX -fastexp 0 global.fa > vf3.tree 
/usr/bin/time -o vf4.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf4.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext AVX -fastexp 2 global.fa > vf4.tree 
/usr/bin/time -o vf5.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf5.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext AVX2 -fastexp 0 global.fa > vf5.tree 
/usr/bin/time -o vf6.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf6.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext AVX2 -fastexp 2 global.fa > vf6.tree 
```

here's what we get for the execution time:

1: 9:01:27
2: 9:04:37
3: 8:50:10
4: 8:44:37
5: 8:46:24
6: 8:44:15

So it's the final settings by a meaningless whisker. Now let's run the analyses

```
/usr/bin/time -o vf1threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf1threads.log -intree iqtree_seqsadded_mp.treefile -threads 1 -double-precision -ext AVX2 -fastexp 2 global.fa > vf1threads.tree 
/usr/bin/time -o vf2threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf2threads.log -intree iqtree_seqsadded_mp.treefile -threads 2 -double-precision -ext AVX2 -fastexp 2 global.fa > vf2threads.tree 
/usr/bin/time -o vf3threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf3threads.log -intree iqtree_seqsadded_mp.treefile -threads 3 -double-precision -ext AVX2 -fastexp 2 global.fa > vf3threads.tree 
/usr/bin/time -o vf4threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf4threads.log -intree iqtree_seqsadded_mp.treefile -threads 4 -double-precision -ext AVX2 -fastexp 2 global.fa > vf4threads.tree 
/usr/bin/time -o vf5threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf5threads.log -intree iqtree_seqsadded_mp.treefile -threads 5 -double-precision -ext AVX2 -fastexp 2 global.fa > vf5threads.tree 
/usr/bin/time -o vf6threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf6threads.log -intree iqtree_seqsadded_mp.treefile -threads 6 -double-precision -ext AVX2 -fastexp 2 global.fa > vf6threads.tree 
/usr/bin/time -o vf7threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf7threads.log -intree iqtree_seqsadded_mp.treefile -threads 7 -double-precision -ext AVX2 -fastexp 2 global.fa > vf7threads.tree 
/usr/bin/time -o vf8threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf8threads.log -intree iqtree_seqsadded_mp.treefile -threads 8 -double-precision -ext AVX2 -fastexp 2 global.fa > vf8threads.tree 
/usr/bin/time -o vf9threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf9threads.log -intree iqtree_seqsadded_mp.treefile -threads 9 -double-precision -ext AVX2 -fastexp 2 global.fa > vf9threads.tree 
/usr/bin/time -o vf10threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf10threads.log -intree iqtree_seqsadded_mp.treefile -threads 10 -double-precision -ext AVX2 -fastexp 2 global.fa > vf10threads.tree 
/usr/bin/time -o vf20threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf20threads.log -intree iqtree_seqsadded_mp.treefile -threads 20 -double-precision -ext AVX2 -fastexp 2 global.fa > vf20threads.tree 
/usr/bin/time -o vf50threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf50threads.log -intree iqtree_seqsadded_mp.treefile -threads 50 -double-precision -ext AVX2 -fastexp 2 global.fa > vf50threads.tree 
/usr/bin/time -o vf100threads.txt -v VeryFastTree -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log vf100threads.log -intree iqtree_seqsadded_mp.treefile -threads 100 -double-precision -ext AVX2 -fastexp 2 global.fa > vf100threads.tree 
```

## Results



