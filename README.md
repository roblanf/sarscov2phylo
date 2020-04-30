# Phylogenetics for SARS-CoV-2

This repository contains scripts for reproducible analyses of SARS-CoV-2 data using GISAID

For many analyses, running your own phylogenies won't be necessary - one can simply use the excellent NextStrain pipeline (particularly augur, for the phylogenetics).

We wanted to build something that did a little more than augur though, which means first solving a few problems in different ways, hence this repository. We hope that if the methods prove useful, we'll have time to build them into augur at a later date.


# Quickstart

1. Get the latest GISAID data (full length and high coverage) in a single fasta file, e.g. `gisaid.fa`

2. Trim the UTRs like this:

`sh trim_seqs.sh -i gisaid.fa -o trimmed.fa -t 8`

3. Align the K most dissimilar sequences:

`sh align_k_dissimilar.sh -i trimmed.fa -o aln_k.fa -k 100`

In the final step we first exclude GISAID sequence with >10 ambiguities, because this helps build a higher-quality alignment of the k most dissimilar sequences. Otherwise we end up selecting the lower quality sequences.


Next steps

* fix profile-alignment of the rest of the seuqences to this alignment (currently semi-working in global_profile_alignment.sh)
* filter the global alignment for bad sites and sequences
* build a global tree