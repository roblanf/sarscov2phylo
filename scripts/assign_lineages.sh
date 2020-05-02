helpFunction()
{
   echo "Assign pangolin lineages to a set of sequences"
   echo "Usage: $0 -i input_fasta -o output_csv -t threads"
   echo "\t-i Full path to fasta file of SARS-CoV-2 sequences (aligned or unaligned)"
   echo "\t-l url to latest lineage file from https://github.com/hCoV-2019/lineages/blob/master/lineages/data/"
   echo "\t-o Output file path for csv file from pangolin"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:l:o:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      i ) lineageurl="$OPTARG" ;;
      o ) outputcsv="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# First we get the pre-computed lineages from hcov lineages

wget $lineageurl > lineages.csv

# get a list of all seqs in inputfasta that are not already in the lineages.csv file
inputdir=$(dirname $inputfasta)
grep '>' $inputfasta  > all_names.txt
cut -d '|' -f2 all_names.txt > all_epiIDs.txt
cut -d',' -f2 lineages.csv > completed_epiIDs.txt
awk 'NR == FNR{ a[$0] = 1;next } !a[$0]' completed_epiIDs.txt all_epiIDs.txt  > epiIDs_todo.txt
grep -F -f epiIDs_todo.txt all_names.txt | tr -d \> > names_todo.txt
faSomeRecords $inputfasta names_todo.txt seqs_todo.fa

# assign lineage with pangolin
pangolin seqs_todo.fa -o $inputdir -t $threads

# TODO combine lineage information into single file
# i.e. add pangolin output to the lineages.csv, being aware they have different columns...
# and noting that lineages.csv is based on epiIDs