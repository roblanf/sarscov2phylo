# Phylogenetics for SARS-CoV-2

This repository contains scripts for reproducible phylogenetic analyses of SARS-CoV-2 data from GISAID.

For many analyses, running your own phylogenies won't be necessary - one can simply use the excellent NextStrain pipeline (particularly augur, for the phylogenetics).

We wanted to build something that did a little more than augur though, which means first solving a few problems in different ways, hence this repository. We hope that if the methods prove useful, we'll have time to build them into augur at a later date.


# Quickstart

1. Get the latest GISAID data (full length and high coverage) in a single fasta file, e.g. `gisaid.fa`, you can do this via the batch download feature on GISAID.

2. Trim the UTRs like this:

`sh trim_seqs.sh -i gisaid.fa -o trimmed.fa -t 8`

where `t` is the number of thread you have

3. Align the K most dissimilar sequences:

`sh align_k_dissimilar.sh -i trimmed.fa -o aln_k.fa -k 100`

In this step we first exclude GISAID sequences with >10 ambiguities, because this helps build a higher-quality alignment of the k most dissimilar sequences. Otherwise we end up selecting the lower quality sequences. Alternatively at this step one would (probably preferably) choose a set of reliable and diverse genome sequences that you know, and that encompass most of the existing diversity in the SARS-CoV-2 data. Regardless, it is highly recommended to check that `aln_k.fa` is a good alignment before proceeding and edit it by hand if it needs any adjustment. The accuracy of all downstream steps depends completely on the accuracy of `aln_k.fa`.

4. Create the global alignment using `aln_k.fa` as a reference alignment 

`sh global_profile_alignment.sh -i trimmed.fa -r aln_k.fa -o global.fa -t 8`

This step aligns each of the sequences in `trimmed.fa` and **ADDS** them to `global.fa`. It's set up like this so that one can easily add new sequences to the global alignment as they are released on GISAID, without re-doing the whole process.  

Next steps

* filter the global alignment for bad sites and sequences 
* build a global tree
* implement global exclusion of sequences via excluded_sequences.txt