#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o phylogenetic_tree -t threads -k num_dissimilar_seqs -a additional_sequences"
   echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -o Output file path for final alignment"
   echo "    -t number of threads to use"
   echo "    -k Number of most dissimilar sequences to align to make the initial guide alignment (suggest ~100)"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:k:a:d:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      k ) k="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$threads" ] || [ -z "$k" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

# first we trim the sequences
inputdir=$(dirname $inputfasta)

cd $inputdir

allseqs=$inputdir"allseqs_unaligned.fasta"
cat $inputfasta > $allseqs


echo ""
echo "Processing raw data and trimming UTRs "
echo ""

trimmed_gisaid="$inputdir/trimmed.fa"
bash $DIR/trim_seqs.sh -i $allseqs -o $trimmed_gisaid -t $threads

#### BUILD THE GLOBAL ALIGNMENT ######

# first align the k most dissimilar sequences, a low k is sensible here
# e.g. 100. Higher numbers slow down the global alignment step
echo ""
echo "Making alignment of $k most dissimilar sequences"
echo ""
aln_k="$inputdir/aln_k.fa"
bash $DIR/align_k_dissimilar.sh -i $trimmed_gisaid -k $k -o $aln_k -t $threads

echo ""
echo "Filtering sites with >10% gaps from k most dissimilar sequence alignment"
echo ""
aln_k_filtered="$inputdir/aln_k_filtered.fa"
esl-alimask --gapthresh 0.1 --informat afa --outformat afa --dna -o $aln_k_filtered -g  $aln_k


echo ""
echo "Making global profile alignment"
echo ""
aln_global="$inputdir/aln_global_unfiltered.fa"
bash $DIR/global_profile_alignment.sh -i $trimmed_gisaid -o $aln_global -t $threads -r $aln_k_filtered


echo ""
echo "Filtering alignment"
echo ""

echo "alignment stats before filtering"
esl-alistat $aln_global

echo ""
echo "Filtering sites with >5% gaps"
esl-alimask --gapthresh 0.05 --informat afa --outformat afa --dna -o $aln_global"_alimask.fa" -g  $aln_global
echo "Filtering sequences that are shorter than 28000 bp and/or have >1000 ambiguities"
esl-alimanip --lmin 28000 --xambig 1000 --informat afa --outformat afa --dna -o $aln_global"_alimanip.fa" $aln_global"_alimask.fa"

mv $aln_global"_alimanip.fa" $outputfasta

echo ""
echo "alignment stats after filtering"
esl-alistat $outputfasta


#### ESTIMATE THE GLOBAL TREE ######

echo ""
echo "Estimating trees with bootstraps using fasttree"
echo ""

# finally, we estimate a tree with 100 bootstraps, using rapidnj, MP, and fasttree
#bash $DIR/tree_nj.sh -i $outputfasta -t $threads
#bash $DIR/tree_mp.sh -i $outputfasta -t $threads
bash $DIR/tree_ft.sh -i $outputfasta -t $threads


echo ""
echo "Cleaning trees with treeshrink"
echo ""

run_treeshrink.py -t $outputfasta'_ft_TBE.tree' -q 0.05 -c -o treeshrink_TBE
run_treeshrink.py -t $outputfasta'_ft_FBP.tree' -q 0.05 -c -o treeshrink_FBP


echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot 'treeshrink_TBE/'$outputfasta'_ft_TBE_0.05.tree' "'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05'" > ft_TBE.tree
nw_reroot 'treeshrink_FBP/'$outputfasta'_ft_FBP_0.05.tree' "'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05'" > ft_FBP.tree

# remove quotes that treeshrink adds
sed -i.bak "s/'//g" ft_TBE.tree
rm ft_TBE.tree.bak

sed -i.bak "s/'//g" ft_FBP.tree
rm ft_FBP.tree.bak

