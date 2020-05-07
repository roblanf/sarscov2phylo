#!/bin/bash

helpFunction()
{
   echo "build a rapidnj tree with 100 bootstraps, support in TBE and FBP"
   echo "Output will be three files:"
   echo "\t *_nj_replicates.tree: 100 bootstrap replicate nj trees"
   echo "\t *_nj_boot_TBE.tree, rapidnj tree with transfer bootstrap supports"
   echo "\t *_nj_boot_FBP.tree, rapidnj tree with felsenstein bootstrap supports"
   echo "Usage: $0 -i fasta_alignment -t threads"
   echo "\t -i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "\t -t number of threads to use"
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

export INPUT_FASTA=$inputfasta


# make the reference tree
echo ""
echo "Making the reference tree with rapidnj"
echo ""
rapidnj $inputfasta -i fa -c $threads -n -t d -x $inputfasta'_rapidnj.tree'


one_bootstrap(){

   bootpre='boot'$1
   goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o $bootpre
   rapidnj $bootpre'0.fa' -i fa -c 1 -n -t d -x $bootpre'.tree'
   sed -i.bak "s/'//g" $bootpre'.tree'
   rm $bootpre'.tree.bak'
   rm $bootpre'0.fa'

}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with rapidnj"
echo ""
boot_nums=($(seq 1 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*.tree > $inputfasta"/_nj_replicates.tree"
inputdir=$(dirname $inputfasta)
find $inputdir -maxdepth 1 -name "boot*" -delete

# make felsenstein bs in iqtre like: iqtree -t TREES_SET_FILE -sup FOCAL_TREE
echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxml-ng --support --tree $inputfasta'_rapidnj.tree' --bs-trees $inputfasta"/_nj_replicates.tree" --prefix $focaltree'nj_boot' --threads $threads --bs-metric fbp,tbe --redo

mv $focaltree'nj_boot.raxml.supportFBP' $inputfasta'_nj_boot_FBP.tree'
mv $focaltree'nj_boot.raxml.supportTBE' $inputfasta'_nj_boot_TBE.tree'

rm $inputfasta'_rapidnj.tree'
rm 'nj_boot.raxml.log'