#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o phylogenetic_tree -t threads -k num_dissimilar_seqs -a additional_sequences"
   echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -o Output file path for final alignment"
   echo "    -t number of threads to use"
   echo "    -k Number of most dissimilar sequences to align to make the initial guide alignment (suggest ~100)"
   echo "    -a full path to unaligned fasta file of additional sequences not on GISAID"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:k:a:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      k ) k="$OPTARG" ;;
      a ) addseqs="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$threads" ] || [ -z "$k" ] || [ -z "$addseqs" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

# first we trim the sequences
inputdir=$(dirname $inputfasta)

cd $inputdir

allseqs=$inputdir"allseqs_unaligned.fasta"
cat $inputfasta $addseqs > $allseqs

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

echo "Filtering sites with >1% gaps"
esl-alimask --gapthresh 0.01 --informat afa --outformat afa --dna -o $aln_global"_alimask.fa" -g  $aln_global
echo "Filtering sequences that are shorter than 29100 bp and/or have >200 ambiguities"
esl-alimanip --lmin 29100 --xambig 200 --informat afa --outformat afa --dna -o $aln_global"alimanip.fa" $aln_global"_alimask.fa"

echo "Ensuring additional seuqences are in the alignment"
grep '>' $addseqs | tr -d \> | faSomeRecords $aln_global"_alimask.fa" /dev/stdin $addseqs"_aln.fa"
cat $aln_global"alimanip.fa" $addseqs"_aln.fa" > $outputfasta"_dupe.fa"
faFilter -uniq $outputfasta"_dupe.fa" $outputfasta

echo "alignment stats after filtering"
esl-alistat $outputfasta


echo ""
echo "Estimating trees with bootstraps"
echo ""

# finally, we estimate a tree with 100 bootstraps, using rapidnj and MP
bash $DIR/tree_nj.sh -i $outputfasta -t $threads
bash $DIR/tree_mp.sh -i $outputfasta -t $threads
bash $DIR/tree_ft.sh -i $outputfasta -t $threads

