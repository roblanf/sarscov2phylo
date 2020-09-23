# Comparing online to de-novo phylogenetics

Since 14-08-20 I have moved to online phylogenetics and I've been logging the fasttree output as I go. The question I want to address here is how online phylogenetics (where I add new sequences every 2 days with MP, then update the resulting tree with ML) compares to standard phylogenetics (where you re-estimate the tree from scratch every time), in terms of the final likelihoods. 

To answer that, all I need to do is re-run every analysis since 14-08-20 with standard phylogenetics. I'll use the following almost default setup in `fasttree` to do that:

```
fasttree -nt -gamma -boot 100 -sprlength 1000 -log fasttree.log global.fa > global.tree
```

This means I can compare execution time (hence wanting the 100 bootstraps to get the SH supports, which matches what I do in the online phylo) and log likelihoods (hence the gamma) in online versus de-novo phylogenetics (nobody calls it de-novo phylogenetics, but it makes sense here to distinguish it from online). 

# Results


