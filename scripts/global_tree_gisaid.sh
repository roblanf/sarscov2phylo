#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o phylogenetic_tree -t threads -k num_dissimilar_seqs -a additional_sequences"
   echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -o Output file path for final alignment"
   echo "    -t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

# first we trim the sequences
inputdir=$(dirname $inputfasta)

cd $inputdir

echo ""
echo "Cleaning raw data"
echo ""

cleaned_gisaid=$inputfasta"_cleaned.fa"
bash $DIR/clean_gisaid.sh -i $inputfasta -o $cleaned_gisaid -t $threads


#### BUILD THE GLOBAL ALIGNMENT ######

echo ""
echo "Making global profile alignment"
echo ""
aln_global="$inputdir/aln_global_unmasked.fa"
bash $DIR/global_profile_alignment.sh -i $cleaned_gisaid -o $aln_global -t $threads


echo ""
echo "Masking alignment"
echo ""
aln_global_masked="$inputdir/aln_global_masked.fa"
bash $DIR/mask_alignment.sh -i $aln_global -o $aln_global_masked -t $threads


echo ""
echo "Filtering sequences that are shorter than 28000 bp and/or have >1000 ambiguities"
echo ""
aln_global_filtered="$inputdir/aln_global_filtered.fa"
esl-alimanip --lmin 28000 --xambig 1000 --informat afa --outformat afa --dna -o $aln_global_filtered $aln_global_masked

echo ""
echo "Removing sites that are >50% gaps, after converting N's to gaps"
echo ""

cp $aln_global_filtered tmp.aln
sed -i.bak '/^[^>]/s/N/-/g' tmp.aln
rm tmp.aln.bak

esl-alimask --gapthresh 0.5 --informat afa --outformat afa --dna -o $outputfasta -g  tmp.aln

rm tmp.aln

echo "alignment stats of global alignment" >> alignments.log
esl-alistat $aln_global >> alignments.log
echo "alignment stats of global alignment after masking sites" >> alignments.log
esl-alistat $aln_global_masked >> alignments.log
echo "alignment stats after filtering out short/ambiguous sequences" >> alignments.log
esl-alistat $aln_global_filtered >> alignments.log
echo "alignment stats of global alignment after trimming sites that are >50% gaps" >> alignments.log
esl-alistat $outputfasta >> alignments.log


#### ESTIMATE THE GLOBAL TREE ######

echo ""
echo "Estimating trees with bootstraps using fasttree"
echo ""

# finally, we estimate a tree with 100 bootstraps using fasttree
bash $DIR/tree_ft.sh -i $outputfasta -t 34


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

echo "Filtering stats"

echo "Unfiltered"
esl-alistat $aln_global
echo "After filtering sites with alimask"
esl-alistat $aln_global"_alimask.fa"
echo "After filtering sequences with alimanip"
esl-alistat $outputfasta
echo "After filtering sequences with TreeShrink"
nw_stats ft_TBE.tree

# zip up for easy file transfer
xz -e -T $threads $outputfasta
xz -e -T $threads $aln_global
xz -e -T $threads $aln_global_filtered
xz -e -T $threads $aln_global_masked
xz -e -T $threads $aln_global_unmasked
xz -e -T $threads $inputfasta
xz -e -T $threads $outputfasta"_ft_replicates.tree"

