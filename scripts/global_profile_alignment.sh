#!/bin/bash

helpFunction()
{
   echo "Align a set of sequences in parallel via global alignment to a reference alignment"
   echo "Usage: $0 -i unaligned_fasta_path -o output_file -r reference_alignment -t threads"
   echo "\t-i Full path to unaligned fasta file of SARS-CoV-2 sequences"
   echo "\t-o Output file path for aligned fasta file. New sequences will be *APPENDED* to this file"
   echo "\t-r Reference alignment (unaligned sequences will be profile-aligned to this)"   
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:o:r:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      r ) refaln="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$refaln" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

set -uxeo pipefail
# so we can get at it as a global variable
export REFERENCE_ALN=$refaln


echo ""
echo "Splitting unaligned input sequences into individual files"
N=$(grep ">" $inputfasta | wc -l)
N=$(( N*2 )) # doubling it is to ensure each record goes to one file
faSplit sequence $inputfasta $N $TMP/profile_align_individual_seq


echo ""
echo "Profile aligning each new sequence to the target alignment"
echo "This can take a while, be patient"
echo ""

profile_align()
{

	seqfile=$1
	alfile=$seqfile"_profile_aligned.fa"
	final=$seqfile"_ind_aligned.fa"

	mafft --thread 1 --quiet --keeplength --add $seqfile "$REFERENCE_ALN" > $alfile

	name=$(grep ">" $seqfile | tr -d ">")

	echo "$name" | faSomeRecords $alfile /dev/stdin $final

	rm $seqfile
	rm $alfile
}

export -f profile_align

find $TMP -type f -name '*profile_align_individual_seq*' | parallel -j $threads --bar "profile_align {}" > /dev/null

# join it all together and clean up
# note that here we *APPEND* to the global alignment, which allows us to add small batches of new stuff whenever we like
find $TMP -name \*.fa_ind_aligned.fa -exec cat {} \; >> "$TMP/aligned_withdups.fa"
find $TMP -maxdepth 1 -name "*profile_align_individual_seq*" -delete

#finally, remove duplicates by name - can occur if a seq. exists in the reference alignment
faFilter -uniq "$TMP/aligned_withdups.fa" $outputfasta
