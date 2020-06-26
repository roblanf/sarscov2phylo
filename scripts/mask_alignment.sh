#!/bin/bash

helpFunction()
{
   echo "Mask a GISAID alignment assuming coordinates are w.r.t. WH1"
   echo "Usage: $0 -i fasta_file_of_SARS-CoV-2_seqs -o output_file -t threads"
   echo "\t-i Full path to input alignment of SARS-CoV-2 sequences"
   echo "\t-o Output file path of masked alignment"
   echo "\t-t number of threads to use"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$threads" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi



# download latest masking file and extract the sites to mask
wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/subset_vcf/problematic_sites_sarsCov2.mask.vcf
mask_sites=$(bcftools query -f '%POS\t' problematic_sites_sarsCov2.mask.vcf)

# TODO: switch this back to just conda when v 0.3.3c is available
# we also need the latest (non-conda) version of goalign for this
wget https://github.com/evolbioinfo/goalign/releases/download/v0.3.3c/goalign_amd64_linux
chmod +x goalign_amd64_linux


# Split GISAID data into individual sequences
# fasplit
echo "Splitting sequences from input file into individual files"
N=$(grep ">" $inputfasta | wc -l)
N=$(( N*2 )) # doubling it is to ensure each record goes to one file
faSplit sequence $inputfasta $N individual_seq


echo "Masking any 7bp window with >=2 uninformative characters, and the first and last 30bp of every sequence"
echo "These sites tend to have a lot more errors than other sites"

export DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Also masking all sites suggested here: https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/subset_vcf/problematic_sites_sarsCov2.mask.vcf"

# build the mask sites list for seqmagick
m=""
for site in $mask_sites; do m="$m,$site:$site"; done
m="${m:1}" # cut the first comma

export MASK_SITES=$m


mask_a_seq()
{

   alfile=$1

   # first we edit in place to mask the first and last 30 characters
   # and any 7bp window with >=2 uninformative characters
   python "$DIR"'/mask_seq.py' $alfile $alfile

   # next we mask the dodgy sites specifically
   seqmagick mogrify --mask $MASK_SITES $alfile 


}

export -f mask_a_seq

# we use parallel like this to avoid "Argument list too long" issues
inputdir=$(dirname $inputfasta)
ls $inputdir | grep individual_seq | parallel -j $threads --bar "mask_a_seq {}" > /dev/null

# make the new sequence set, and clean up
# use find to avoid "Argument list too long" issue
echo "putting all the trimmed sequences into one file"
find $inputdir -name \*individual_seq*.fa -exec cat {} \; > $outputfasta

#clean up
echo "cleaning up"
find $inputdir -maxdepth 1 -name "individual_seq*" -delete