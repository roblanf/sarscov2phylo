#!/bin/bash

helpFunction()
{
   echo "Align the K most dissimilar sequences from a file"
   echo "Usage: $0 -i unaligned_fasta_path -k num_seqs -o output_file -t threads"
   echo "\t-i Full path to unaligned fasta file of alignable sequences"
   echo "\t-k Number of most dissimilar sequences to align"
   echo "\t-o Output file path"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:k:o:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      k ) k="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$k" ] || [ -z "$outputfasta" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

n=$(grep '>' $inputfasta | wc -l)
echo ""
echo ""
echo "Filtering the $n initial unaligned input sequences to remove any with > 50 ambiguous bases"
echo ""
echo ""
seqmagick quality-filter --max-ambiguous 50 --max-length 100000 --min-length 100 $inputfasta $inputfasta"_filtered.fa"
inputfasta=$inputfasta"_filtered.fa"

# how many sequences in the inputfasta
n=$(grep '>' $inputfasta | wc -l)
echo "$n sequences remain after filtering"

#define a cutoff to switch to the quick & dirty method of choosing the seqs
cutoff=50000

inputdir=$(dirname $inputfasta)

echo "Using MAFFT to get a guide tree"

if (( $n > $cutoff )); then

	echo "Using parttree distance method in MAFFT"

	# use this method for maximum speed on v large datasets
	mafft --thread $threads --retree 0 --treeout --parttree $inputfasta > /dev/null

	# replace all newlines
	tr -d '\n' < "$inputfasta.tree" > "$inputdir/guide1.tree"

	# add branchlengths because this method doesn't give any.
	# assuming equal branchlengths is obviously not correct, but nevertheless
	# still gives a pragmatic way of selecting the k divergent sequences
	sed 's/),/)0.1,/g' "$inputdir/guide1.tree" > "$inputdir/guide2.tree"
	sed 's/))/)0.1)/g' "$inputdir/guide2.tree" > "$inputdir/guide3.tree"

else

	echo "Using 6mer distance method in MAFFT"

	# use this method if it's feasible, which calculates 6mer distances
	mafft --thread $threads --retree 0 --treeout $inputfasta > /dev/null
	# replace all newlines
	tr -d '\n' < $inputfasta.tree > "$inputdir/guide3.tree"

fi

# now continue for all methods 

# add a semicolon to the tree so iq-tree can read it
echo ";" >> "$inputdir/guide3.tree"

echo "Using IQ-TREE to select the K most dissimilar sequences"

# get the k most dissimilar sequences using iq-tree
iqtree -te "$inputdir/guide3.tree" -k $k

echo "Extracting K most dissimilar sequences"

# this part just gets the list of names that we want, and keeps it around
# as k_names.txt
if (( $n > $cutoff )); then
	# extract those from the fasta file to make a new fasta
	# https://www.biostars.org/p/2822/
	# get the numbers of the sequences in the sorted file
	awk '/^The optimal PD set has/{flag=1;next}/^Corresponding sub-tree/{flag=0}flag' "$inputdir/guide3.tree.pda" > "$inputdir/numbers.txt"
	# get the names
	grep '>' $inputfasta > "$inputdir/all_names.txt"
	# get the names given the line numbers
	awk 'FNR==NR{a[$1];next}(FNR in a){print}' "$inputdir/numbers.txt" "$inputdir/all_names.txt" > "$inputdir/k_names_raw.txt"
	# remove the ">"
	tr -d \> < "$inputdir/k_names_raw.txt" > "$inputdir/k_names.txt"

else
	# extract those from the fasta file to make a new fasta
	awk '/^The optimal PD set has/{flag=1;next}/^Corresponding sub-tree/{flag=0}flag' "$inputdir/guide3.tree.pda" > "$inputdir/k_names_long.txt"

	# these names are like 123_name. For safety we only use the number, since mafft edits special characters in names
	cut -d _ -f 1 "$inputdir/k_names_long.txt" > "$inputdir/numbers.txt"

	# get the names
	grep '>' $inputfasta > "$inputdir/all_names.txt"

	# get the names given the line numbers
	awk 'FNR==NR{a[$1];next}(FNR in a){print}' "$inputdir/numbers.txt" "$inputdir/all_names.txt" > "$inputdir/k_names_raw.txt"

	# remove the ">"
	tr -d \> < "$inputdir/k_names_raw.txt" > "$inputdir/k_names.txt"

fi


rm "$inputdir/k_names_long.txt"
rm "$inputdir/guide1.tree"
rm "$inputdir/guide2.tree"
rm "$inputdir/guide3.tree"
rm "$inputdir/guide3.tree.log"
rm "$inputfasta.tree"
rm "$inputdir/k_names_raw.txt"
rm "$inputdir/numbers.txt"
rm "$inputdir/all_names.txt"
rm "$inputdir/guide3.tree.pda"


# extract the fasta records we want
#seqtk subseq $inputfasta "$inputdir/k_names.txt" > "$inputdir/k_unaligned.fa"
faSomeRecords $inputfasta "$inputdir/k_names.txt" "$inputdir/k_unaligned.fa"

echo "Aligning K most dissimilar sequences"

# align the k most dissimilar sequences in MAFFT
mafft --thread $threads --auto "$inputdir/k_unaligned.fa" > $outputfasta

rm "$inputdir/k_unaligned.fa"
rm $inputfasta
rm "$inputdir/k_names.txt"
