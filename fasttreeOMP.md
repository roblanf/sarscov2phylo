# Benchmarking fasttree OMP

Question - how much faster is fasttree OMP (compiled on my machine) than the non-OMP version I get from conda, and how many threads is optimal.

This one's simple. First I get fasttree OMP and compile it:

```
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm
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


## Results



