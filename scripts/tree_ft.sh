#!/bin/bash

helpFunction()
{
   echo "build a fasttree tree with 100 bootstraps, support in TBE and FBP"
   echo "Usage: $0 [-o outprefix] -i fasta_alignment -t threads"
   echo "    -i Full path to aligned fasta file of SARS-CoV-2 sequences"
   echo "    -o Prefix for output trees"
   echo "    -t number of threads to use"
   echo "Output will be three files:"
   echo "    *_ft_replicates.tree: 100 bootstrap replicate fasttree trees"
   echo "    *_ft_boot_TBE.tree, fasttree tree with transfer bootstrap supports"
   echo "    *_ft_boot_FBP.tree, fasttree tree with felsenstein bootstrap supports"
   exit 1 # Exit script after printing help
}

while getopts "i:t:o:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outprefix="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
outprefix="${outprefix:-${inputfasta%%.fa}}"
if [ ! -d $outprefix ]
then
	# add an underscore to output prefix if it's not a directory
	outprefix="${outprefix}_"
fi

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$threads" ] || [ -z "${outprefix}" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

set -uxeo pipefail

export TMP="${TMP:-$(mktemp -d)}"
export INPUT_FASTA=$inputfasta
DIR="$(cd "$(dirname "$0")" && pwd)"
base="$(basename "$inputfasta" .fa)"

# make the reference tree
echo ""
echo "Making the reference tree with fasttree -fastest option"
echo ""
fasttree -nosupport -nt -fastest $inputfasta > "$TMP/${base}_multi.fasttree"

Rscript $DIR/bifurcate.R -i "$TMP/${base}_multi.fasttree" -o "$TMP/${base}.fasttree"

one_bootstrap(){
	base="$(basename "$INPUT_FASTA" .fa)"
	bootpre="$TMP/boot_${base}_$1"
	goalign build seqboot -i "$INPUT_FASTA" -t 1 -n 1 -S -o "$bootpre"
	fasttree -nosupport -nt -fastest "${bootpre}.fa" > "${bootpre}_unrooted.tree"
	nw_reroot "${bootpre}_unrooted.tree" 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > "${bootpre}_multi.tree"
}

export -f one_bootstrap

echo ""
echo "Making 100 bootstrap trees with fasttree -fastest"
echo ""
boot_nums=($(seq 1 100))
parallel -j $threads --bar "one_bootstrap {}" ::: ${boot_nums[@]} > /dev/null

# make the file we need and clean up
find "$TMP" -maxdepth 1 -name "boot_${base}*multi.tree" | xargs cat > "$TMP/${base}_ft_replicates_multi.tree"

# make it bifurcating
Rscript $DIR/bifurcate.R -i "$TMP/${base}_ft_replicates_multi.tree" -o "$TMP/${base}_ft_replicates.tree" 

find $TMP -maxdepth 1 -name "boot_${base}_*" -delete

echo ""
echo "Running raxml to map bootstrap support to focal tree"
raxpre="$TMP/${base}_ft_boot"
raxml-ng --support --tree "$TMP/${base}.fasttree" --bs-trees "$TMP/${base}_ft_replicates.tree" --prefix "$raxpre" --threads $threads --bs-metric fbp,tbe --redo

# re-root as per https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1
echo ""
echo "Re-rooting tree on hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
echo "see https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1"
echo ""
nw_reroot "${raxpre}.raxml.supportFBP" 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > "${outprefix}ft_FBP.tree"
rm "${raxpre}.raxml.supportFBP"

nw_reroot "${raxpre}.raxml.supportTBE" 'hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05' > "${outprefix}ft_TBE.tree"
rm "${raxpre}.raxml.supportTBE"

rm "${raxpre}.raxml.log"
