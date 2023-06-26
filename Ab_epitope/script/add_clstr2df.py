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

def add_clstr2df(cluster_identity):
    # add cluster into data table

    with open('result/cluster.tsv','w') as o:
        o.write('\t'.join(['Name',f'cluster_ID_by_{c}'])+'\n')
        for record in SeqIO.parse(f"result/cluster/epitope_all{c}.fasta.clstr", "fasta"):
            c_id = record.id
            s = record.seq
            uniprot_ids = [x.split('...')[0] for x in s.split('>') if '...' in x]
            for i in uniprot_ids:
                o.write('\t'.join([str(i),c_id]) + '\n')
    # add original file so we can add clstr to this file

    metadata_DF = pd.read_excel('result/epitope_info_clstr_v2.xlsx')
    cluster_df = pd.read_csv('result/cluster.tsv', sep='\t')
    DF=pd.merge(metadata_DF,cluster_df,on='Name')
    DF.to_excel('result/epitope_info_clstr_v2.xlsx',index=False)
    os.remove('result/cluster.tsv')

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='Add cluster to an Excel table')
    parser.add_argument('-i','--input_file', default='data/epitope_info_v2.xlsx', type=str, help='the input Excel file')
    parser.add_argument('-o','--output_file',  default='result/epitope_info_clstr_v2.xlsx',type=str, help='the output xlsx file')
    parser.add_argument('-l', '--list',  default=['0.8','0.9','0.95','0.99','0.995'], nargs='+', help='cluster identity list to add the xlsx table')
    # Parse the arguments
    args = parser.parse_args()
    shutil.copy(args.input_file, args.output_file)
    for c in args.list:
        preprocess_clstr(f"result/cluster/epitope_all{c}.fasta.clstr")
        add_clstr2df(cluster_identity=c)


