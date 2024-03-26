import pandas as pd
import argparse

if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='Add cluster to an Excel table')
    parser.add_argument('-i','--input_OASfile', default='data/raw_data/OAS_memory_paired.csv', type=str, help='the input csv file')
    parser.add_argument('-o','--output_file',  default='mBLM/result/memory_paired_Abs.csv',type=str, help='the output csv file')
    
    parser.add_argument('-gb','--genbank', default='data/raw_data/all_paired_antibodies_from_GB_v6.xlsx',type=str, help='the paired antibody from genbank')
    
    # Parse the arguments
    args = parser.parse_args()
    OAS_DF = pd.read_csv(args.input_OASfile, low_memory=False)
    genbank_df = pd.read_excel(args.genbank)
    # merge main table and the cluster id by the chosen cluster
    # Concatenate the dataframes based on the same columns ['B', 'C']
    concatenated = pd.concat([OAS_DF, genbank_df])
    deduplicated_df = concatenated.drop_duplicates(subset=['sequence_alignment_aa_heavy','sequence_alignment_aa_light'], keep='first')
    deduplicated_df = deduplicated_df.drop_duplicates(subset=['Name'], keep='first')
    # filter heavy chain length by 100
    filtered_df = deduplicated_df[deduplicated_df['sequence_alignment_aa_heavy'].str.len() > 100]
    filtered_df.to_csv(args.output_file,index=False)
