#!/bin/bash

helpFunction()
{
   echo "Clean sequences from GISAID according to known issues"
   echo "Usage: $0 -i fasta_file_of_SARS-CoV-2_seqs -o output_file -t threads"
   echo "\t-i Full path to unaligned fasta file of SARS-CoV-2 sequences"
   echo "\t-o Output file path"
   echo "\t-t number of threads to use"
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


input_seqs=$inputfasta
output_seqs=$outputfasta

echo ""
echo "Preprocessing GISAID data to fix known issues"

# first let's fix cases where gisaid doesn't put a newline before the '>''
echo "Making sure all records start with a new line before the '>'"
sed 's/>/\'$'\n>/g' $input_seqs > $input_seqs"_namesfixed.fa"
input_seqs=$input_seqs"_namesfixed.fa"

sed -i.bak '/^$/d' $input_seqs


# next we remove sequences in the excluded_sequences.tsv file
echo "Removing sequences in exluded_sequence.tsv"
BASEDIR=$(dirname "$0")
exseq=$BASEDIR/../excluded_sequences.tsv
cut -f1 $exseq | faSomeRecords $input_seqs /dev/stdin $input_seqs"tmp.fa" -exclude
rm $input_seqs
mv $input_seqs"tmp.fa" $input_seqs

# now we fix spaces in the fasta headers from GISAID
# we do this ONLY on header lines
echo "Replacing spaces in sequence names with '_'"
sed -i.bak '/^>/s/ /_/g' $input_seqs

echo "Replacing brackets in sequence names with '-'"
sed -i.bak '/^>/s/(/-/g' $input_seqs
sed -i.bak '/^>/s/)/-/g' $input_seqs

# now we replaces spaces in the sequences with N's
# these occur in a number of PHE UK sequences
echo "Replacing spaces in sequences with 'N'"
sed -i.bak '/^[^>]/s/ /N/g' $input_seqs
rm $input_seqs'.bak'

echo ""


mv $input_seqs $outputfasta

