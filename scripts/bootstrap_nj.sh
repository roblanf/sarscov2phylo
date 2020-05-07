#!/bin/bash

helpFunction()
{
   echo "Bootstrap 100 rapidnj analyses using IQ-TREE"
   echo "Usage: $0 -i fasta_alignment -o output_file -t threads"
   echo "\t-i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "\t-f Full path focal phylogeny on which to map bootstrap support"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:f:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      f ) focaltree="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$focaltree" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

export INPUT_FASTA=$inputfasta


one_bootstrap(){

   bootpre='boot'$1
   goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o $bootpre
   rapidnj $bootpre'0.fa' -i fa -c 1 -n -t d -x $bootpre'.tree'
   sed -i.bak "s/'//g" $bootpre'.tree'
   rm $bootpre'.tree.bak'
   rm $bootpre'0.fa'

}

export -f one_bootstrap

boot_nums=($(seq 1 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*.tree > nj_replicates.tree
inputdir=$(dirname $inputfasta)
find $inputdir -maxdepth 1 -name "boot*" -delete

# make felsenstein bs in iqtre like: iqtree -t TREES_SET_FILE -sup FOCAL_TREE
echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxml-ng --support --tree $focaltree --bs-trees nj_replicates.tree --prefix $focaltree'nj_boot' --threads $threads --bs-metric fbp,tbe --redo
