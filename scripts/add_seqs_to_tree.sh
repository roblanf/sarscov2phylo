#!/bin/bash

helpFunction()
{
   echo "Input a tree and an alignment - remove sequences from the tree that are not in the alignment, then add sequences to the tree that are in the alignment but not the tree"
   echo "Usage: $0 -i fasta_file_of_SARS-CoV-2_seqs -o output_file -t threads"
   echo "\t-i Full path to input alignment of SARS-CoV-2 sequences"
   echo "\t-s Full path to input starting tree of SARS-CoV-2 sequences"
   echo "\t-o Full path to output tree, which will contain all sequences in the alignment"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:s:o:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      i ) inputtree="$OPTARG" ;;
      o ) outputtree="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$inputtree" ] || [ -z "$outputtree" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


# 1. Get the sequence names from the alignment

grep ">" $inputfasta | cut -c 2- > alignment_names.txt

# 2. Remove sequence from the tree that are not in the alignment
# this can happen e.g. if a sequence was removed from GISAID, or excluded because it was on a long branch
# or if the filtering conditions changed

