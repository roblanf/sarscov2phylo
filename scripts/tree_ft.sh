#!/bin/bash

helpFunction()
{
   echo "build a fasttree tree with 100 bootstraps, support in TBE and FBP"
   echo "Usage: $0 -i fasta_alignment -t threads"
   echo "    -i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "    -t number of threads to use"
   echo "Output will be three files:"
   echo "    *_ft_replicates.tree: 100 bootstrap replicate fasttree trees"
   echo "    *_ft_boot_TBE.tree, fasttree tree with transfer bootstrap supports"
   echo "    *_ft_boot_FBP.tree, fasttree tree with felsenstein bootstrap supports"
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
DIR="$(cd "$(dirname "$0")" && pwd)"


one_bootstrap(){

   if (( $(bc -l <<< "$1 == 0") )); then
      # make the reference tree if bootstrap ref = 0
      # parameters are set to be more thorough for the ML tree than for the bootstrap trees
      # this is pure pragmatism. For reasoning behind parameters, see tree_estimation2.md in this repo
      echo ""
      echo "Making the reference tree with fasttree -fastest option"
      echo ""
      fasttree -nosupport -nt -gamma -sprlength 20 -refresh 0.8 -topm 1.5 -close 0.75 "$INPUT_FASTA" > "$INPUT_FASTA"'multi.fasttree'

   else

      bootpre='boot'$1
      goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o $bootpre
      fasttree -nosupport -nt -fastest $bootpre'0.fa' > $bootpre'unrooted.tree'

   fi   
}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with fasttree -fastest"
echo ""
boot_nums=($(seq 0 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*unrooted.tree > $inputfasta"_ft_replicates_multi.tree"

inputdir=$(dirname $inputfasta)

find $inputdir -maxdepth 1 -name "boot*" -delete

# we hardcode 55 threads here, to keep it fast
gotree compute support tbe -i $inputfasta"multi.fasttree" -b $inputfasta"_ft_replicates_multi.tree" -t 55 -o $inputfasta"_ft_TBE_unrooted.tree"
gotree compute support fbp -i $inputfasta"multi.fasttree" -b $inputfasta"_ft_replicates_multi.tree" -t 55 -o $inputfasta"_ft_FBP_unrooted.tree"

# re-root as per https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1
echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot $inputfasta"_ft_FBP_unrooted.tree" 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ft_FBP.tree'
rm $inputfasta"_ft_FBP_unrooted.tree"

nw_reroot $inputfasta"_ft_TBE_unrooted.tree" 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ft_TBE.tree'
rm $inputfasta"_ft_TBE_unrooted.tree"

