
helpFunction()
{
   echo "Filter strains and columns from an alignment"
   echo "Usage: $0 -i fasta_alignment -o output_fasta"
   echo "    -i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "    -o Output file path for filtered fasta"   
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
echo "Filtering strains shorter than 29100 and/or with more than 400 ambiguities (that's ~2%)"
echo "Filtering sites with more than 1% gaps"
echo ""

esl-alimask --gapthresh 0.01 --informat afa --outformat afa --dna -o $inputfasta"_alimask.fa" -g  $inputfasta
esl-alimanip --lmin 29100 --xambig 400 --informat afa --outformat afa --dna -o $outputfasta $inputfasta"_alimask.fa"

rm $inputfasta"_alimask.fa"