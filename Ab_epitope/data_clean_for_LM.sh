#!/bin/bash
python mBLM/script/merge_raw_data.py
python mBLM/script/df2fasta.py mBLM/result/memory_paired_Abs.csv mBLM/result/memory_paired_Abs.fasta

# clustering
output_prefix=mBLM/result/cluster/
f=mBLM/result/memory_paired_Abs.fasta
name=memory_paired_Abs
for x in 0.5 0.6
do
	cd-hit -i $f -o $output_prefix$name$x -c $x -M 32000 -d 0 -T 32 -n 3 -aL 0.8 -s 0.95  -uS 0.2  -sc 1 -sf 1
done
    
for x in 0.7 0.8 0.9 0.95
do
	cd-hit -i $f -o $output_prefix$name$x -c $x -M 32000 -d 0 -T 32 -n 5 -aL 0.8 -s 0.95  -uS 0.2  -sc 1 -sf 1
done  


# add cluster id to dataset

python mBLM/script/add_clstr2df.py -i mBLM/result/memory_paired_Abs.csv -o mBLM/result/memory_paired_Abs_final.csv

# split dataset
python mBLM/script/split_dataset.py -i mBLM/result/memory_paired_Abs_final.csv -o data/dataset/memory_paired_Abs
