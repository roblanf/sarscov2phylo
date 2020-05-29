#!/bin/bash

helpFunction()
{
   echo "Use ML to refine the sub-tree of a focal taxon"
   echo "Usage: $0 -i global_alignment -t threads -g global_phylogeny -f focal_taxon_name -d depth"
   echo "    -i Full path to global aligned fasta file of SARS-CoV-2 sequences"
   echo "    -t number of threads to use"
   echo "    -g global phylogeny that contains the focal taxon"
   echo "    -f focal taxon name (must be present in alignment and global phylogeny)"
   echo "    -d number of substitutions below the focal taxon to prune the sub-tree"
   exit 1 # Exit script after printing help
}

while getopts "i:t:g:f:d:" opt
do
   case "$opt" in
      i ) aln="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      g ) tree="$OPTARG" ;;
      f ) seq="$OPTARG" ;;
      d ) depth="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$aln" ] || [ -z "$threads" ] || [ -z "$tree" ] || [ -z "$seq" ] || [ -z "$depth" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


echo ""
echo "Calculating number of nodes to go back from focal sequence to cut tree"
echo ""

# calculate the minimum depth we need to get to 
length=$(esl-alistat $aln | grep 'Alignment length' | cut -d ':' -f2 | tr -d '[:space:]')
depth_bl=$( bc -l <<< "$depth/$length")


# get the depth we need to go to in order to get the right number of substitutions
# we just keep going one parent deeper until the distance from the focal leaf
# to the parent node is at least depth/alignment_length
# we have two additional conditions
# number of taxa must be >100
# break if number of taxa >999 (because we can't pracically estimate trees that large...)
f=0
for c in $(seq 1 100); do

    depth=$(nw_clade -c $c $tree $seq | nw_distance - $seq)

    ntax=$(nw_clade -c $c $tree $seq | nw_labels -I - | wc -l)

    # convert depth to a float
    printf -v depth "%f" "$depth"

    if (( $(bc -l <<< "$depth > $depth_bl") )); then
        if (( $(bc -l <<< "$ntax > 100") )); then
            break
        fi
    fi

    if (( $(bc -l <<< "$ntax > 1000") )); then
        
        # keep datasets to a maximum of 1000 taxa
        c=$(( $c - 1 ))
        f=1 # just flag that we did this for later...
        break

    fi

done

echo ""
echo "Cutting tree "$c" nodes deeper than the focal sequence"
echo ""


echo ""
echo "Making a sub-alignment"
echo ""
# make a sub-alignment
# the nw_ commands just get the sub-tree, get the names
# then we cut and pass to faSomeRecords
nw_clade -c $c $tree $seq | nw_distance -n -s f - | cut -f1 | faSomeRecords $aln /dev/stdin $seq'_aln.fa'

# add the outgroup sequence to the alignment if it's not already there

og="hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05"
ogn=$(grep $og $seq'_aln.fa' | wc -l)
if (( $ogn == 0 )); then

    echo $og | faSomeRecords $aln /dev/stdin $aln'WH4.fa'
    cat $seq'_aln.fa' $aln'WH4.fa' > tmp.fa
    mv tmp.fa $seq'_aln.fa'
    rm tmp.fa
    rm $aln'WH4.fa'

fi


echo ""
echo "Sub-alignment statistics"
echo ""
esl-alistat $seq'_aln.fa'


# run sub-alignment to get tree.
DIR="$(cd "$(dirname "$0")" && pwd)"

# get the ML tree with bootstraps
bash $DIR/tree_ft.sh -i $seq'_aln.fa' -t $threads

finalTBE=$seq'_aln.fa_ft_TBE.tree'
finalFBP=$seq'_aln.fa_ft_FBP.tree'

# make rudimentary pdfs
echo ""
echo "Attempting to build PDFs of your trees (if inkscape works)"
echo ""


echo "fill:blue L "$seq > css.map
echo '"stroke-width:2; stroke:blue"  Clade '$seq >> css.map

nw_display -s -w 1000 -c css.map -i 'font-size:8; fill:red' -b 'font-size:8; fill:green' $finalTBE > $finalTBE'.svg'
inkscape -f $finalTBE'.svg' -D -A $finalTBE.pdf

nw_display -s -w 1000 -c css.map -i 'font-size:8; fill:red' -b 'font-size:8; fill:green' $finalFBP > $finalFBP'.svg'
inkscape -f $finalFBP'.svg' -D -A $finalFBP.pdf

rm css.map

# dating the tree
#../iqtree-2.0.4-Linux/bin/iqtree2 -s global.fa -te global.fa_mp_boot_TBE.tree -keep-ident -m JC -fixbr --date TAXNAME -o "hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05" -nt 4 -pre date

