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

   if (( $(bc -l <<< "$e == 0") )); then
      # make the reference tree if bootstrap ref = 0
      echo ""
      echo "Making the reference tree with fasttree -fastest option"
      echo ""
      fasttree -nosupport -nt -fastest "$INPUT_FASTA" > "$INPUT_FASTA"'multi.fasttree'

      Rscript $DIR/bifurcate.R -i "$INPUT_FASTA"'multi.fasttree' -o "$INPUT_FASTA"'.fasttree'

   else

      bootpre='boot'$1
      goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o $bootpre
      fasttree -nosupport -nt -fastest $bootpre'0.fa' > $bootpre'unrooted.tree'
      nw_reroot $bootpre'unrooted.tree' "'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05'" > $bootpre'multi.tree'

   fi   
}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with fasttree -fastest"
echo ""
boot_nums=($(seq 0 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*multi.tree > $inputfasta"_ft_replicates_multi.tree"

# make it bifurcating
Rscript $DIR/bifurcate.R -i $inputfasta"_ft_replicates_multi.tree" -o $inputfasta"_ft_replicates.tree" 

inputdir=$(dirname $inputfasta)

find $inputdir -maxdepth 1 -name "boot*" -delete
rm $inputfasta"_ft_replicates_multi.tree"

echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxml-ng --support --tree $inputfasta'.fasttree' --bs-trees $inputfasta"_ft_replicates.tree" --prefix $inputfasta'_ft_boot' --threads $threads --bs-metric fbp,tbe --redo

# re-root as per https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1
echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot $inputfasta'_ft_boot.raxml.supportFBP' 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ft_FBP.tree'
rm $inputfasta'_ft_boot.raxml.supportFBP'

nw_reroot $inputfasta'_ft_boot.raxml.supportTBE' 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ft_TBE.tree'
rm $inputfasta'_ft_boot.raxml.supportTBE'

rm $inputfasta'_ft_boot.raxml.log'
