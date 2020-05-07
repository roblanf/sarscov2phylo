#!/bin/bash

helpFunction()
{
   echo "build a stepwise addition parsimony tree with 100 bootstraps, support in TBE and FBP"
   echo "Usage: $0 -i fasta_alignment -t threads"
   echo "Output will be three files:"
   echo "\t *_mp_replicates.tree: 100 bootstrap replicate mp trees"
   echo "\t *_mp_boot_TBE.tree, stepwise mp tree with transfer bootstrap supports"
   echo "\t *_mp_boot_FBP.tree, stepwise mp tree with felsenstein bootstrap supports"
   echo "\t-i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

inputdir=$(dirname $inputfasta)
export INPUT_FASTA=$inputfasta

# make the reference tree
echo ""
echo "Making the reference tree with IQ-TREE stepwise addition parsimony"
echo ""
iqtree -s $inputfasta -keep-ident -n 0 -m JC -fixbr -nt 1 



one_bootstrap(){

   bootpre='boot'$1
   iqtree -s "$INPUT_FASTA" -keep-ident -n 0 -m JC -fixbr -nt 1 -bo 1 -wbtl -pre $bootpre

}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with IQ-TREE stepwise addition parsimony"
echo ""
boot_nums=($(seq 1 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*.treefile > $inputfasta"mp_replicates.tree"
inputdir=$(dirname $inputfasta)
find $inputdir -maxdepth 1 -name "boot*" -delete

# make felsenstein bs in iqtre like: iqtree -t TREES_SET_FILE -sup FOCAL_TREE
echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxml-ng --support --tree $inputfasta'.parstree' --bs-trees $inputfasta"_mp_replicates.tree" --prefix $focaltree'mp_boot' --threads $threads --bs-metric fbp,tbe --redo

mv $focaltree'mp_boot.raxml.supportFBP' $inputfasta'_mp_boot_FBP.tree'
mv $focaltree'mp_boot.raxml.supportTBE' $inputfasta'_mp_boot_TBE.tree'

rm $inputfasta.*
rm $inputdir'/mp_boot.raxml.log'