import argparse
import pandas as pd

# usage: python df2fasta.py input_excel output.fasta
if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='Convert an Excel table to a FASTA file.')
    parser.add_argument('input_file', metavar='input_file', type=str, help='the input Excel file')
    parser.add_argument('output_file', metavar='output_file', type=str, help='the output FASTA file')
    parser.add_argument('-name', default='Name',type=str, help="Name column")
    parser.add_argument('-seq', default='sequence_alignment_aa_heavy',type=str, help="sequence column")
    # Parse the arguments
    args = parser.parse_args()

    # Read the Excel file into a Pandas DataFrame
    df = pd.read_csv(args.input_file)

    # Open the output file for writing
    with open(args.output_file, 'w') as f:
        # Loop over the rows in the DataFrame
        for index, row in df.iterrows():
        # Write the header line with the sequence name
            f.write('>{}\n'.format(row[args.name]))
        # Write the sequence to the file
            f.write('{}\n'.format(row[args.seq]))

