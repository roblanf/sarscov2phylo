# new approach

# the script shoudl be almost the same as the global phylogeny script

# with the sole difference that a starting tree is given

# this tree should contain at least some of the same sequence names as given in the alignment

# the process here is fairly simple (I think...)

# 1. Use R to trim the input tree down to only sequences found in the final aligned fasta file (global.fa)

# get names from alignment


# 2. Use IQ-TREE with the tree from (1) as a starting tree, using MP to place the new sequences

# 3. Use the tree from 2 to input as a starting tree to fasttree


# Of course, this makes bootstraps redundant. So there are two options here:

# 1. Just use fasttree SH supports (seems fine for now)
# 2. Try IQ-TREE aLRT supports (issue: uses a lot of memory, possibly slow? let's see)

 