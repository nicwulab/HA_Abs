from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import shutil
import argparse
import os

'''
ADD CLUSTER ID INTO MAIN TABLE
'''

def preprocess_clstr(clstr_file):
    '''
    modify cluster file to be readable as fasta
    :param clstr_file:
    :return: auto replace old file
    '''

    # Open the cluster file for reading
    with open(clstr_file, "r") as f:
        # Initialize empty lists to store the sequence headers and sequences
        headers = []
        sequences = []
        # Initialize empty variables to store the current header and sequence
        current_header = ""
        current_sequence = ""
        # Loop through each line in the file
        for line in f:
            # If the line starts with ">", it's a sequence header
            if line.startswith(">"):
                # Replace any spaces in the header with underscores
                current_header = line.strip().replace(" ", "_")
                headers.append(current_header)
                # If there was a previous sequence, add it to the sequences list
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            # If the line doesn't start with ">", it's part of the current sequence
            else:
                current_sequence += line
        # Add the final sequence to the sequences list
        sequences.append(current_sequence)

    # write out the updated cluster file with the same name
    with open(clstr_file, "w") as f:
        # Loop through the headers and sequences lists and write them to the output file
        for i in range(len(headers)):
            f.write(headers[i] + "\n")
            f.write(sequences[i])

def add_clstr2df(tmp,output_path,cluster):
    # add cluster into data table

    with open(f'{tmp}/cluster.tsv','w') as o:
        o.write('\t'.join(['Name',f'cluster_ID_by_{cluster}'])+'\n')
        for record in SeqIO.parse(f"mBLM/result/cluster/memory_paired_Abs{cluster}.clstr", "fasta"):
            c_id = record.id
            s = record.seq
            # loop over all the sequences id in the same cluster, find the representative using '... *'
            uniprot_ids = [x.split('...')[0] for x in s.split('>') if '...' in x]
            for i in uniprot_ids:
                o.write('\t'.join([str(i),c_id]) + '\n')
    # add original file so we can add clstr to this file

    metadata_DF = pd.read_csv(output_path, low_memory=False)
    cluster_df = pd.read_csv(f'{tmp}/cluster.tsv', sep='\t')
    # merge main table and the cluster id by the chosen cluster
    DF=pd.merge(metadata_DF,cluster_df,on='Name', how='left')
    DF.to_csv(output_path,index=False)
    os.remove(f'{tmp}/cluster.tsv')

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='Add cluster to an Excel table')
    parser.add_argument('-i','--input_file', default='mBLM/result/memory_paired_Abs.csv', type=str, help='the input csv file')
    parser.add_argument('-o','--output_file',  default='mBLM/result/memory_paired_Abs_final.csv',type=str, help='the output csv file')
    parser.add_argument('-l', '--list',  default=['0.5','0.6','0.7','0.8','0.9','0.95'], nargs='+', help='cluster identity list to add the xlsx table')
    parser.add_argument('-tmp','--tmp_dir', default='mBLM/result',type=str, help='the output cluster file temporary directory')
    
    # Parse the arguments
    args = parser.parse_args()
    shutil.copy(args.input_file, args.output_file)
    for c in args.list:
        preprocess_clstr(f"mBLM/result/cluster/memory_paired_Abs{c}.clstr")
        add_clstr2df(args.tmp_dir,args.output_file,c)


