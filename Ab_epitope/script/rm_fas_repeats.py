import argparse
from Bio import SeqIO
def remove_duplicates(input_file, output_file):
    # Create a dictionary to store the sequences
    seen_ids = set()
    seq_dict = {}

    # Open the input and output files
    with open(input_file, "r") as f_in:
        # Loop over the sequences in the input file
        for record in SeqIO.parse(f_in, "fasta"):
            # Check if the sequence identifier has been seen before
            if str(record.id) not in seen_ids:
                # If not, add it to the set of seen identifiers and write the sequence to the output file
                seen_ids.add(str(record.id))
                seq_dict[str(record.id)]=str(record.seq)
    with open(output_file, "w") as f_out:
        for name,seq in seq_dict.items():

            f_out.write(f'>{name}\n')
            f_out.write(f'{seq}\n')



if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='Remove repeat name sequence from a FASTA file.')
    parser.add_argument('input_file', metavar='input_file', type=str, help='the input Excel file')
    parser.add_argument('output_file', metavar='output_file', type=str, help='the output FASTA file')
    # Parse the arguments
    args = parser.parse_args()

    remove_duplicates(args.input_file,args.output_file)
