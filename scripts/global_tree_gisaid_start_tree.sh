#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o output_filename -s start_tree -t threads "
   echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -p Filepath to previous iteration to be updated (must contain an ft_SH.tree file)"
   echo "    -t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:p:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      p ) previous="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$previous" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

inputdir=$(dirname $inputfasta)
outputfasta=global.fa
inputtree=previous_iteration_files/ft_SH.tree # this is copied over in a few lines

cd $inputdir

# first we copy over the previous iterations files for reproducibility
mkdir previous_iteration_files
cp $previous'/ft_SH.tree' previous_iteration_files/

# first we trim the sequences, but we don't remove any
echo ""
echo "Cleaning raw data"
echo ""

cleaned_gisaid=$inputfasta"_cleaned.fa"
bash $DIR/clean_gisaid.sh -i $inputfasta -o $cleaned_gisaid -t $threads


#### BUILD THE GLOBAL ALIGNMENT ######

echo ""
echo "Making global profile alignment"
echo ""
aln_global="$inputdir/1_aln_global_unmasked.fa"
bash $DIR/global_profile_alignment.sh -i $cleaned_gisaid -o $aln_global -t $threads


echo ""
echo "Masking alignment"
echo ""
aln_global_masked="$inputdir/2_aln_global_masked.fa"
bash $DIR/mask_alignment.sh -i $aln_global -o $aln_global_masked -t $threads


echo ""
echo "Filtering sequences that are shorter than 28000 bp and/or have >1000 ambiguities"
echo ""
aln_global_filtered="$inputdir/3_aln_global_filtered.fa"
esl-alimanip --lmin 28000 --xambig 1000 --informat afa --outformat afa --dna -o $aln_global_filtered $aln_global_masked

echo ""
echo "Removing sites that are >50% gaps, after converting N's to gaps"
echo ""
cp $aln_global_filtered tmp.aln
sed -i.bak '/^[^>]/s/N/-/g' tmp.aln
rm tmp.aln.bak
aln_4="$inputdir/4_aln_global_gapmasked.fa"
esl-alimask --gapthresh 0.5 --informat afa --outformat afa --dna -o $aln_4 -g  tmp.aln

rm tmp.aln

echo "sequences downloaded from GISAID" >> alignments.log
echo $(grep ">" $inputfasta | wc -l) >> alignments.log
echo "//" >> alignments.log
echo "alignment stats of global alignment" >> alignments.log
esl-alistat $aln_global >> alignments.log
echo "alignment stats of global alignment after masking sites" >> alignments.log
esl-alistat $aln_global_masked >> alignments.log
echo "alignment stats after filtering out short/ambiguous sequences" >> alignments.log
esl-alistat $aln_global_filtered >> alignments.log
echo "alignment stats of global alignment after trimming sites that are >50% gaps" >> alignments.log
esl-alistat $aln_4 >> alignments.log


#### ESTIMATE THE GLOBAL TREE ######


# 1. Take the alignment and the tree, remove the irrelevant sequences with R, then add the new seuqences with IQ-TREE

echo ""
echo "Removing unused sequences from input tree"
echo ""
grep ">" $aln_4 | cut -c 2- > alignment_names.txt
Rscript $DIR/clean_tree.R $inputtree alignment_names.txt

echo ""
echo "Adding new sequences to input tree with IQ-TREE"
echo ""
# get the latest IQ-TREE
wget https://github.com/iqtree/iqtree2/releases/download/v2.1.0/iqtree-2.1.0-Linux.tar.gz
tar -xvzf iqtree-2.1.0-Linux.tar.gz

# this just adds the new sequences with parsimony
# benchmarking shows that 1 thread is optimal
./iqtree-2.1.0-Linux/bin/iqtree2 -seed 1729 -s $aln_4 -g input_tree_cleaned.tree -n 0 -m JC -fixbr -nt 1 --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -pre iqtree_seqsadded_mp

echo "Tree stats after adding sequences with IQ-TREE" >> alignments.log
nw_stats iqtree_seqsadded_mp.treefile >> alignments.log


# 2. Clean the tree by removing long terminals then running treeshrink

# Iteratively remove terminal branches that are >=5 mutations long until there are none such
# this also writes de-novo an exlucded_sequences.tsv file
echo ""
echo "Cleaning tree by removing long terminal branches"
echo ""
# memory insurance, was hitting odd issues on CentOS that seemed to be due to this
ulimit -s 64000
Rscript $DIR/prune_terminals.R iqtree_seqsadded_mp.treefile iqtree_pruned_terminals.tree

echo "Tree stats after pruning branches with >5 mutations" >> alignments.log
nw_stats iqtree_pruned_terminals.tree >> alignments.log


# Run TreeShrink to identify sequences on long branches that were not already pruned
echo ""
echo "Cleaning tree with treeshrink"
echo ""
run_treeshrink.py -t iqtree_pruned_terminals.tree -q 0.05 -c -o treeshrink

echo "Tree stats after pruning branches with >5 mutations" >> alignments.log
nw_stats treeshrink/iqtree_pruned_terminals_0.05.tree >> alignments.log

# now we update the excluded sequences file
echo ""
echo "Updating excluded sequences file with sequences exluded by TreeShrink"
echo ""
Rscript $DIR/update_excluded_seqs.R excluded_sequences.tsv treeshrink_SH/iqtree_pruned_terminals_RS_0.05.txt



exseq=excluded_sequences.tsv
cut -f1 $exseq | faSomeRecords $aln_4 /dev/stdin $outputfasta -exclude

echo "alignment stats of global alignment after removing sequences with treeshrink" >> alignments.log
esl-alistat $outputfasta >> alignments.log


# Run FastTree on the new alignmnet and the tree from TreeShrink...
echo ""
echo "Optimising treeshrink-ed tree with fasttree MP"
echo ""

# we have to do some contortions to set the optimal number of threads for fasttree, which is 3 (see fasttreeOMP.md)
env > old_env.txt
old_threads=$(grep -hoP '^OMP_NUM_THREADS=\K\d+' old_env.txt)
rm old_env.txt
export OMP_NUM_THREADS=3

FastTreeMP -nt -gamma -nni 0 -spr 2 -sprlength 1000 -boot 100 -log fasttree.log -intree treeshrink/iqtree_pruned_terminals_0.05.tree $outputfasta > $outputfasta'_ft_SH.tree'

if [ -n "$old_threads" ]; then
    export OMP_NUM_THREADS=$old_threads
else
    unset OMP_NUM_THREADS
fi

echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot $outputfasta'_ft_SH.tree' "'EPI_ISL_406801'" > ft_SH.tree

sed -i.bak "s/'//g" ft_SH.tree
rm ft_SH.tree.bak

echo "Tree stats after re-optimising tree with FastTree" >> alignments.log
nw_stats ft_SH.tree >> alignments.log

echo "//"
echo "Number of new sequences added this iteration" >> alignments.log
wc -l alignment_names_new.txt >> alignments.log


# zip up for easy file transfer
xz -e -T $threads $aln_4
xz -e -T $threads $outputfasta
xz -e -T $threads $aln_global
xz -e -T $threads $aln_global_filtered
xz -e -T $threads $aln_global_masked
xz -e -T $threads $aln_global_unmasked
xz -e -T $threads $inputfasta
xz -e -T $threads $inputfasta"_cleaned.fa"
xz -e -T $threads fasttree.log
xz -e -T $threads iqtree_seqsadded_mp.iqtree

rm goalign_amd64_linux
rm -rf iqtree-2.1.0-Linux/
rm iqtree_seqsadded_mp.uniqueseq.phy
