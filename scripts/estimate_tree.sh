
helpFunction()
{
   echo "Estimate a global tree for a large alignment of SARS-CoV-2 sequences"
   echo "Usage: $0 -i fasta_alignment -o output_tree -t threads"
   echo "\t-i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "\t-o Output file path for phylogenetic tree"
   echo "\t-b number of bootstrap trees you want"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:" opt
do 
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) output_tree="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$output_tree" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


echo ""
echo "Estimating a tree for your data"
echo ""


echo "Estimating a MP tree in IQ-TREE"
iqtree -s $inputfasta -keep-ident -n 0 -m JC -fixbr -nt 1 -pre pars
iqtree -s $inputfasta -keep-ident -n 0 -m JC -fixbr -nt 7 -b0 7 -wbt -pre parsboot

# get the focal tree

# parallel to 
   # generate a bootstrap alignment from the focal alignment
   # run it through rapidnj/IQ-TREE on one core
   # write tree file

# cat the tree files together (check for newlines)

# iq-tree to map the bootstrap trees to the focal tree