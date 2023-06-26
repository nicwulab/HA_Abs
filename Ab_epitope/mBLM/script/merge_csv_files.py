import os
import pandas as pd
import argparse

def merge_csv_files(input_dir, output_file, delimiter):

    # Create the empyt list of data frames
    dfs = []

    # Loop through all csv.gz files in the directory
    for filename in os.listdir(input_dir):
        if filename.endswith('.csv.gz'):

            # Read the csv.gz file into a pandas DataFrame
            df = pd.read_csv(os.path.join(input_dir, filename), header=1, compression='gzip', delimiter=delimiter)
            # List to hold the column names
            column_names = ['sequence_id_heavy', 'sequence_heavy','sequence_alignment_aa_heavy', 'v_call_heavy','d_call_heavy','j_call_heavy', 'cdr1_aa_heavy','cdr2_aa_heavy','cdr3_aa_heavy','sequence_id_light','sequence_light','sequence_alignment_aa_light', 'v_call_light','d_call_light','j_call_light', 'cdr1_aa_light','cdr2_aa_light','cdr3_aa_light']
            df = df[column_names]
            
            # Add the DataFrame to the list
            dfs.append(df)
            #merge dataframe
            merge_df=pd.concat(dfs)

    # Write the merged DataFrame to a CSV file
    merge_df['Name'] = merge_df['sequence_id_heavy'] + '_' + merge_df['v_call_heavy']
    merge_df.to_csv(output_file, index=False)

    print("Merged CSV file has been created!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge a set of csv.gz files into a single csv file.')
    parser.add_argument('input_dir', help='Path to the input directory containing the csv.gz files.')
    parser.add_argument('output_file', help='Name of the output csv file.')
    parser.add_argument('--delimiter', default=',', help='Column delimiter for the csv files. Default is ",".')

    args = parser.parse_args()

    merge_csv_files(args.input_dir, args.output_file, args.delimiter)

