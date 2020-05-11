


# input the taxa you want to refine a tree for

# also the global trees you want to use

# also the depth you want to go from each focal taxon

# get the tree N levels deep - could just keep increasing N until we have at least 100 sequences

nw_clade -c 5 global.fa_mp_boot_TBE.tree 'ACT01_consensus_sequence'

# to get how many taxa 
nw_clade -c 50 global.fa_mp_boot_TBE.tree 'ACT01_consensus_sequence' | nw_stats -

# to get how deep we've gone from the focal leaf
nw_clade -c 40 global.fa_mp_boot_TBE.tree 'ACT01_consensus_sequence' | nw_distance - 'ACT01_consensus_sequence'



depth=3
tree=global.fa.parstree
aln=global.fa
seq=ACT01_consensus_sequence
threads=50

# calculate the minimum depth we need to get to 
length=$(esl-alistat $aln | grep 'Alignment length' | cut -d ':' -f2 | tr -d '[:space:]')
depth_bl=$( bc -l <<< '$depth/$length')


# get the depth we need to go to in order to get the right number of substitutions
# we just keep going one parent deeper until the distance from the focal leaf
# to the parent node is at least depth/alignment_lenght
for c in $(seq 1 100); do
	depth=$(nw_clade -c $c $tree $seq | nw_distance - $seq)

	# convert depth to a float
	printf -v depth "%f" "$depth"

	if (( $(bc -l <<< "$depth > $depth_bl") )); then
		break
	fi
done


# make a sub-alignment
# the nw_ commands just get the sub-tree, get the names
# then we cut and pass to faSomeRecords

nw_clade -c $c $tree $seq | nw_distance -n -s f - | cut -f1 | faSomeRecords $aln /dev/stdin $seq'_aln.fa'

# add the outgroup sequencde to the alignment
echo "hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05" | faSomeRecords $aln /dev/stdin $aln'WH4.fa'

cat $seq'_aln.fa' $aln'WH4.fa' > tmp.fa
mv tmp.fa $seq'_aln.fa'
rm tmp.fa

# run sub-alignment in IQ-TREE with rigorous settings.
DIR="$(cd "$(dirname "$0")" && pwd)"

# get the ML tree with bootstraps
bash $DIR/tree_ml.sh -i $seq'_aln.fa' -t $threads

finalTBE=$seq'_aln.fa_ml_boot_TBE.tree'
finalFBP=$seq'_aln.fa_ml_boot_FBP.tree'

# display the tree
echo '"<circle style=’fill:red;stroke:black’ r=’5’/>" I "'$seq'"' > ornament.map
nw_display -s -o ornament.map $finalTBE > $finalTBE'.svg'
inkscape -f $finalTBE'.svg' -D -A $finalTBE.pdf


# dating the tree
#../iqtree-2.0.4-Linux/bin/iqtree2 -s global.fa -te global.fa_mp_boot_TBE.tree -keep-ident -m JC -fixbr --date TAXNAME -o "hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05" -nt 4 -pre date