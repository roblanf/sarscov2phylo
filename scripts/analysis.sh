#!/bin/bash

helpFunction()
{
   echo "Make an ML phylogeny from a large FASTA file of GISAID sequences"
   echo "Usage: $0 -i GISAID_fasta -o phylogenetic_tree -t threads"
   echo "\t-i Full path to unaligned fasta file of SARS-CoV-2 sequences from GISAID"
   echo "\t-o Output file path for phylogenetic tree"
   echo "\t-t number of threads to use"
   echo "\t-k Number of most dissimilar sequences to align to make the initial constraint tree"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:k:" opt
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

trimmed_gisaid="$inputdir/trimmed_gisaid.fa"
sh $DIR/trim_seqs.sh -i $inputfasta -o $trimmed_gisaid -t $threads

#### BUILD THE GLOBAL ALIGNMENT ######

# first align the k most dissimilar sequences, a low k is sensible here
# e.g. 100. Higher numbers slow down the global alignment step
aln_k="$inputdir/k_aligned.fa"
sh $DIR/align_k_dissimilar.sh -i $trimmed_gisaid -k $k -o $aln_k

global_aln="$inputdir/global_alignment.fa"
sh $DIR/global_profile_alignment.sh -i $trimmed_gisaid -o $global_alignment -t 8 -r $aln_k





##### CALCULATE A TREE #########

# get a bionj tree of all the data

# use it as a starting tree for IQ-TREE model selection

# now update the tree and brlens using that model



# now do a full IQ-TREE run on the K most dissimilar sequences
# based on many previous analyses, we look at GTR models with a range of different rate parameterisations
#iqtree -s $aln_k -m MFP -nt $threads -merit AICc -mset "GTR" -mrate "E,I,G,I+G,R,I+R" -pre k_aligned -bb 10000

# get the best model, we'll use this later
#iqfile="k_aligned.iqtree"
#bestmod=$(grep "Best-fit model according to" $iqfile | cut -d':' -f2 | cut -d' ' -f2)



#### ADD THE REST OF THE SEQUENCES TO THE TREE ####


