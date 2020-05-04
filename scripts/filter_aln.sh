
helpFunction()
{
   echo "Filter strains and columns from an alignment"
   echo "Usage: $0 -i fasta_alignment -o output_fasta -t threads"
   echo "\t-i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "\t-o Output file path for filtered fasta"
   exit 1 # Exit script after printing help
}

while getopts "i:o:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

echo ""
echo "Filtering strains shorter than 29400 and/or with more than 20 ambiguities"
echo "Filtering sites with more than 0.1% gaps"
echo ""

esl-alimanip --lmin 29400 --xambig 20 --informat afa --outformat afa --dna -o $inputfasta"_alimanip.fa" $inputfasta 
esl-alimask --gapthresh 0.001 --informat afa --outformat afa --dna -o $outputfasta -g  $inputfasta"_alimanip.fa"

rm $inputfasta"_alimanip.fa"