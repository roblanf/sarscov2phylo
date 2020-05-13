#!/bin/bash

helpFunction()
{
   echo "build a ML tree in raxml-ng with 100 bootstraps, support in TBE and FBP, uses WH04 as an ougroup in all cases"
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

# make the reference tree
echo ""
echo "Making the reference tree with IQTREE with GTR+I+G model"
echo ""
raxml-ng --msa $inputfasta --model GTR+I+G --threads 4  --blmin 0.0000001 --tree 'pars{5}'

one_bootstrap(){

   bootpre='boot'$1
   goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o $bootpre
   # these are bootstraps, and MP is very good, so we'll just do 1 starting tree
   raxml-ng --msa $bootpre'0.fa' --model GTR+I+G --threads 1 --blmin 0.0000001 --tree 'pars{1}'

}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with raxml-ng with GTR+I+G model"
echo ""
boot_nums=($(seq 1 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
cat boot*.bestTree > $inputfasta"_ml_replicates.tree"
inputdir=$(dirname $inputfasta)
find $inputdir -maxdepth 1 -name "boot*" -delete


echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxml-ng --support --tree $inputfasta'.raxml.bestTree' --bs-trees $inputfasta"_ml_replicates.tree" --prefix $inputfasta'_ml_boot' --threads $threads --bs-metric fbp,tbe --redo

# re-root as per https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1
echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot $inputfasta'_ml_boot.raxml.supportFBP' 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ml_FBP.tree'
rm $inputfasta'_ml_boot.raxml.supportFBP'

nw_reroot $inputfasta'_ml_boot.raxml.supportTBE' 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > $inputfasta'_ml_TBE.tree'
rm $inputfasta'_ml_boot.raxml.supportTBE'

rm $inputfasta'_ml_boot.raxml.log'

rm $inputfasta'.raxml'*
