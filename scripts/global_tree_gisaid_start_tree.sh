#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o output_filename -s start_tree -t threads "
   echo "    -i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -o Output file path for final alignment"
   echo "    -s Full path to input starting tree of SARS-CoV-2 sequences"
   echo "    -t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:o:s:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      s ) inputtree="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$inputtree" ] || [ -z "$outputfasta" ] || [ -z "$threads" ]
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
esl-alistat $outputfasta >> alignments.log


#### ESTIMATE THE GLOBAL TREE ######


# 1. Take the alignment and the tree, remove the irrelevant sequences with R, then add the new seuqences with IQ-TREE

echo ""
echo "Removing unused sequences from input tree"
echo ""
grep ">" $outputfasta | cut -c 2- > alignment_names.txt
Rscript $DIR/clean_tree.R $inputtree alignment_names.txt

echo ""
echo "Adding new sequences to input tree with IQ-TREE"
echo ""
# get the latest IQ-TREE
wget https://github.com/iqtree/iqtree2/releases/download/v2.1.0/iqtree-2.1.0-Linux.tar.gz
tar -xvzf iqtree-2.1.0-Linux.tar.gz

# this just adds the new sequences with parsimony
# benchmarking shows that 1 thread is optimal
./iqtree-2.1.0-Linux/bin/iqtree2 -s $outputfasta -g input_tree_cleaned.tree -n 0 -m JC -fixbr -nt 1 --suppress-zero-distance --suppress-list-of-sequences --suppress-duplicate-sequence -pre iqtree_seqsadded_mp

echo ""
echo "Optimising tree with fasttree"
echo ""
fasttree -nt -gamma -nni 0 -spr 10 -sprlength 1000 -boot 100 -intree iqtree_seqsadded_mp.treefile $outputfasta > $outputfasta'_ft_SH.tree'


echo ""
echo "Cleaning tree with treeshrink"
echo ""
run_treeshrink.py -t $outputfasta'_ft_SH.tree' -q 0.05 -c -o treeshrink_SH


echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot 'treeshrink_SH/'$outputfasta'_ft_SH_0.05.tree' "'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05'" > ft_SH.tree


sed -i.bak "s/'//g" ft_SH.tree
rm ft_SH.tree.bak


echo "After filtering sequences with TreeShrink" >> alignments.log
nw_stats ft_TBE.tree >> alignments.log


# zip up for easy file transfer
xz -e -T $threads $outputfasta
xz -e -T $threads $aln_global
xz -e -T $threads $aln_global_filtered
xz -e -T $threads $aln_global_masked
xz -e -T $threads $aln_global_unmasked
xz -e -T $threads $inputfasta
xz -e -T $threads $inputfasta"_cleaned.fa"

rm goalign_amd64_linux
rm -rf iqtree-2.1.0-Linux/

# tar up the files for easy transfer
tar -zcvf dat.tar.gz *