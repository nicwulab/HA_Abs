#!/bin/bash

# export fasta file
python ./script/df2fasta.py ./data/epitope_info_v2.xlsx ./data/epitope_all.fasta

# clustering

fasta_suffix=.fasta
fasta_prefix=./data/
output_prefix=./result/cluster/
f=epitope_all.fasta
name=`basename $f ".fasta"`
for x in 0.8 0.9 0.92 0.95 0.98 0.99 0.995
do
	cd-hit -i $f -o $output_prefix$name$x$fasta_suffix -c $x -M 32000 -d 0 -T 8 -n 5 -aL 0.8 -s 0.95  -uS 0.2  -sc 1 -sf 1
done


# add cluster id to dataset

python ./script/add_clstr2df.py

# split dataset
python ./script/split_dataset.py 
# remove repeat fasta sequence
pytho script/rm_fas_repeats.py result/epitope_clean.fasta result/epitope_clean_v2.fasta 