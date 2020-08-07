#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o phylogenetic_tree -t threads -k num_dissimilar_seqs -a additional_sequences"
   echo "    -i Full unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "    -p Phylogeny of SARS-CoV-2 sequences passed via -i"
   echo "    -a unaligned fasta file of additional sequences"
   echo "    -t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:p:a:t:" opt
do
   case "$opt" in
      i ) input_aln="$OPTARG" ;;
      p ) input_tree="$OPTARG" ;;
      a ) input_seqs="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$input_aln" ] || [ -z "$input_tree" ] || [ -z "$input_seqs" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

DIR="$(cd "$(dirname "$0")" && pwd)"

# first we trim the input sequences
inputdir=$(dirname $input_aln)

cd $inputdir

echo ""
echo "Cleaning unaligned raw data"
echo ""

cleaned_input_seqs=$input_seqs"_cleaned.fa"
bash $DIR/clean_gisaid.sh -i $input_seqs -o $cleaned_input_seqs -t $threads


#### BUILD THE GLOBAL ALIGNMENT ######

echo ""
echo "Making new sequence profile alignment"
echo ""
aln_seqs_unmasked="$inputdir/aln_seqs_unmasked.fa"
bash $DIR/global_profile_alignment.sh -i $cleaned_input_seqs -o $aln_seqs_unmasked -t $threads


echo ""
echo "Masking new sequences"
echo ""
aln_seqs_masked="$inputdir/aln_seqs_masked.fa"
bash $DIR/mask_alignment.sh -i $aln_seqs_unmasked -o $aln_seqs_masked -t $threads




echo ""
echo "Removing masked sites from new alignment"
echo "Note this may not always be sufficient for adding the new sequences to the global alignment"
echo "In this case, you will need to use mafft or similar to profile-align the two alignment after this step"
echo ""

aln_seqs_removed="$inputdir/aln_seqs_removed.fa"
esl-alimask --gapthresh 0.99 --informat afa --outformat afa --dna -o $aln_seqs_removed -g  $aln_seqs_masked


esl-alistat $aln_seqs_removed
esl-alistat $input_aln

echo ""
echo "Creating total alignment by adding new sequences to global alignment"
echo ""

cat $input_aln $aln_seqs_removed > total.fa


#### ESTIMATE THE GLOBAL TREE ######

echo ""
echo "Estimating position of new sequences with IQ-TREE"
echo ""

iqtree -s total.fa -g ft_TBE.tree -czb -nt 4 --keep-ident -m GTR+G -pers 0.1 -nbest 1 -ninit 1



echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot total.fa.parstree 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > total.fa.pars.treenw_reroot 'treeshrink_FBP/'$outputfasta'_ft_FBP_0.05.tree' "'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05'" > ft_FBP.tree



# zip up for easy file transfer
xz -e -T $threads $input_aln
xz -e -T $threads $input_seqs
xz -e -T $threads $aln_seqs_unmasked
xz -e -T $threads $aln_seqs_removed
xz -e -T $threads $aln_seqs_masked
xz -e -T total.fa


rm goalign_amd64_linux

# tar up the files for easy transfer
tar -zcvf dat.tar.gz *
