#!/bin/bash
fasta_suffix=.fasta
fasta_prefix=./data/
output_prefix=./result/cluster/
for f in ./data/*.fasta
do
	name=`basename $f ".fasta"`
	for x in 0.8 0.9 0.92 0.95 0.98 0.99 0.995
	do
		cd-hit -i $f -o $output_prefix$name$x$fasta_suffix -c $x -M 32000 -d 0 -T 8 -n 5 -aL 0.8 -s 0.95  -uS 0.2  -sc 1 -sf 1
	done
done


