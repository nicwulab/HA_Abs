from transformers import AutoTokenizer, AutoModel
import torch
from Bio import SeqIO
import argparse
import torch.distributed as dist
from tqdm import tqdm
from pathlib import Path

def extract_features(model, tokenizer, fasta_file, output_dir,device, batch_size = 32):
    sequences = []
    fasta_ids = []
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        fasta_id = record.id
        sequences.append(sequence)
        fasta_ids.append(fasta_id)
    
    # Tokenize the sequences
    encoded_inputs = tokenizer(sequences, padding=True, truncation=True, return_tensors="pt")
    # Move the model to the specified device
    model = model.to(device)
    model.eval()
    # Batch processing
    num_sequences = len(sequences)
    num_batches = (num_sequences + batch_size - 1) // batch_size
    
    for batch_idx in tqdm(range(num_batches), desc="Extracting features", leave=False):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, num_sequences)
        
        batch_inputs = {k: v[start_idx:end_idx].to(device) for k, v in encoded_inputs.items()}
        fasta_ids_batch = fasta_ids[start_idx:end_idx]
        
        # Pass the batch through the model
        with torch.no_grad():
            model_outputs = model(**batch_inputs, output_hidden_states=True)
        
        # Extract the desired feature layer, make sure to change the device to cpu
        feature_layer = model_outputs.hidden_states[-1].to(device="cpu")  # Assuming the last layer
        
        # Save individual PT files
        for idx, fasta_id in enumerate(fasta_ids_batch):
            sequence_features = feature_layer[idx]
            output_file = f"{output_dir}/{fasta_id}.pt"
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            torch.save(sequence_features, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract per-token representations and model outputs for sequences in a FASTA file")
    parser.add_argument(
        "--model_location",
        default='./mBLM',
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",)
    parser.add_argument(
        "--fasta_file",
        default='./data/OAS_memory_paired_clean.fasta',
        type=str,
        help="FASTA file on which to extract representations",)
    parser.add_argument(
        "--output_dir",
        default='./data/dataset/mBLM_embedding',
        type=str,
        help="output directory for extracted representations",)
    parser.add_argument("--batch_size", type=int, default=16, help="batch size")
    args = parser.parse_args()
    tokenizer = AutoTokenizer.from_pretrained(args.model_location)
    model = AutoModel.from_pretrained(args.model_location)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    extract_features(model, tokenizer, args.fasta_file, args.output_dir,device,args.batch_size)