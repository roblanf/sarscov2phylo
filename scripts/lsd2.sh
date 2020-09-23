# date with LSD2

git clone https://github.com/tothuhien/lsd2.git
cd lsd2/src/
make

cd ../..

# get the dates
nw_labels -I ft_SH.tree > labels.txt
cut -d '|' -f 3  labels.txt > dates.txt
paste labels.txt dates.txt > datefile.txt

# add number of lines to start of datefile.txt
l=$(wc -l datefile.txt | cut -d' ' -f1)

echo $l | cat - datefile.txt > temp && mv temp datefile.txt

# lsd2
# -i input file
# -s sequence length
# -l 0.000000005 because this is the minimum brlen in the ft_SH.tree
./lsd2/bin/lsd2_unix -i ft_SH.tree -s 29000 -l 0.000000005 -d datefile.txt -v 1 -a 'b(2019-10-30,2019-12-17)'