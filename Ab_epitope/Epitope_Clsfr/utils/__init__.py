import torch
from torch.utils.data import Dataset,DataLoader
import subprocess

from torch.nn import functional as F
from sklearn.preprocessing import MinMaxScaler,LabelEncoder
from multiprocessing import Pool
from functools import partial
import numpy as np
import pandas as pd
import os
import glob
import gzip
import sys
from Bio import SeqIO
import esm


def load_select_terms_and_weight(distribution_file):
    df = pd.read_csv(distribution_file, sep='\t')
    df['weight'] = 1 - (df['Number of pdb'] /
                        df['Number of pdb'].sum())
    select_terms = df.Term.to_list()
    weights = df.weight.to_list()
    return select_terms, torch.FloatTensor(weights)


def seq2onehot(seq):
    """Create 26-dim embedding"""
    chars = ['-', 'D', 'G', 'U', 'L', 'N', 'T', 'K', 'H', 'Y', 'W', 'C', 'P',
             'V', 'S', 'O', 'I', 'E', 'F', 'X', 'Q', 'A', 'B', 'Z', 'R', 'M']
    vocab_size = len(chars)
    vocab_embed = dict(zip(chars, range(vocab_size)))

    # Convert vocab to one-hot
    vocab_one_hot = np.zeros((vocab_size, vocab_size), int)
    for _, val in vocab_embed.items():
        vocab_one_hot[val, val] = 1

    embed_x = [vocab_embed[v] for v in seq]
    seqs_x = np.array([vocab_one_hot[j, :] for j in embed_x])
    return seqs_x


def zero_padding(input_tensor, max_h=None, max_w=None):
    h, w = input_tensor.shape
    h_pad = max_h
    w_pad = max_w
    # create a new tensor of the desired size
    pad_tensor = torch.zeros(h_pad, w_pad)

    # copy the values from the input tensor to the pad tensor
    h, w = min(h_pad, input_tensor.shape[0]), min(w_pad, input_tensor.shape[1])
    pad_tensor[0:h, 0:w] = input_tensor[0:h, 0:w]

    return pad_tensor


def load_npz(npz_path):
    '''
    load trRosetta structural features
    :return:data
    '''
    uni_id = npz_path.split('/')[-1].split(".")[0]
    # print(npz_path)

    npz_dict = {k: v for k, v in np.load(npz_path, allow_pickle=True).items()}

    embedded_seq = npz_dict['embedded_seq']
    # onehot embedding
    onehot_seq = seq2onehot(npz_dict['seqres'].item())
    C_alpha = torch.from_numpy(npz_dict['C_alpha'])
    A = torch.le(C_alpha, 10).long()  # Adjacency matrix
    edge_index = A.nonzero().t().contiguous()
    input_phi12 = npz_dict['trRosetta'].item()['CA1-CB1-CB2']
    input_phi21 = npz_dict['trRosetta'].item()['CA2-CB2-CB1']
    input_theta12 = npz_dict['trRosetta'].item()['N1-CA1-CB1-CB2']
    input_theta21 = npz_dict['trRosetta'].item()['N2-CA2-CB2-CB1']
    input_omega = npz_dict['trRosetta'].item()['CA1-CB1-CB2-CA2']
    y = npz_dict['labels']
    data = Data(x=torch.from_numpy(onehot_seq).float(),
                edge_index=edge_index,
                embedded_seq=torch.from_numpy(embedded_seq),
                input_phi12=zero_padding(torch.from_numpy(input_phi12).nan_to_num(), max_w=1500),
                input_phi21=zero_padding(torch.from_numpy(input_phi21).nan_to_num(), max_w=1500),
                input_theta12=zero_padding(torch.from_numpy(input_theta12).nan_to_num(), max_w=1500),
                input_theta21=zero_padding(torch.from_numpy(input_theta21).nan_to_num(), max_w=1500),
                input_omega=zero_padding(torch.from_numpy(input_omega).nan_to_num(), max_w=1500),
                seq=npz_dict["seqres"].item(),
                ID = npz_dict['id'].item(),
                y=torch.IntTensor([y]))
    return data
def load_esm_embedding(path):

    esm_embedding = torch.load(path)['representations'][33]
    return esm_embedding
def load_mBLM_embedding(path):

    mBLM_embedding = torch.load(path)
    return mBLM_embedding

class EpitopeDataset(Dataset):
    '''
    Custom dataset
    
    self.df is input dataframe with col [Name,VH_AA,Antigen_epitopes...]
    self.embed_path is the path for embedded sequence files(*.pt)
    self.prediction is True only when you are using predict.py
    '''
    def __init__(self, df,embed_path,prediction=False):
        super().__init__()
        self.df = df
        self.embed_path = embed_path
        self.prediction = prediction
        self.label_map={"HA:Head":0, "HA:Stem":1,"HIV":2, "S:NTD":3, "S:RBD":4, "S:S2":5, "Others":6}
        self.load_esm_embedding = partial(load_esm_embedding)
        self.load_mBLM_embedding = partial(load_mBLM_embedding)
    def __len__(self):
         return len(self.df)
    def __getitem__(self, index):
        row = self.df.iloc[index]
        name = row['Name']
        seq = row['VH_AA']
        target_data = row['Antigen_epitopes']
        # Convert your input and target data into PyTorch tensors
        if 'esm2' in self.embed_path:
            esm_embedded_tensor = self.load_esm_embedding(f'{self.embed_path}{name}.pt')
            x=zero_padding(esm_embedded_tensor,max_h=150,max_w=1280)
        elif "mBLM" in self.embed_path:
            embedded_tensor = self.load_mBLM_embedding(f'{self.embed_path}{name}.pt')
            x=zero_padding(embedded_tensor,max_h=150,max_w=768)
        else:
            onehot_seq = seq2onehot(seq)
            x = torch.from_numpy(onehot_seq).float()
            x=zero_padding(x,max_h=150,max_w=26)
        if self.prediction:
            y = np.array(0)
        else:
            y = np.array(self.label_map[target_data])
        return x,y


def df2fasta(df,output_file,name='Name',seq='VH_AA'):
    # Open the output file for writing
    with open(output_file, 'w') as f:
        for index, row in df.iterrows():

            f.write(f'>{row[name]}\n')
            f.write(f'{row[seq]}\n')

def get_dataset_from_df(df_path, batch_size, ncpu=16):
    if df_path.split('.')[-1] == 'tsv':
        df = pd.read_csv(df_path,sep='\t')
    elif df_path.split('.')[-1] == 'csv':
        df = pd.read_csv(df_path)
    elif df_path.split('.')[-1] == 'xlsx':
        df = pd.read_excel(df_path)
    else:
        print('ERROR: the format of dataframe is not defined')
        
    filename = os.path.basename(df_path).split('.')[0]
    # first, the get the fasta file and remove the depulicates
    fas_f = f'result/{filename}.fasta'
    if not os.path.exists(fas_f):
        print('start generating fasta file')
        df2fasta(df,fas_f)
        # remove duplicate sequence from fasta
        subprocess.run(['python', 'script/rm_fas_repeats.py', fas_f, fas_f])
    # run esm embedding
    embedding_output_path = f'result/{filename}_mBLM_embedding/'
    os.makedirs(embedding_output_path, exist_ok=True)
    if len(os.listdir(embedding_output_path)) < 1:
        print('start embedding seqeuence...')
        subprocess.run(['python', 'extract_mBLM_feature.py',  '--model_location','mBLM','--fasta_file',fas_f,'--output_dir', embedding_output_path])
    test_loader = DataLoader(EpitopeDataset(df,embedding_output_path,prediction=True), batch_size=batch_size,shuffle=False,
                                         num_workers=ncpu, pin_memory=True)
    return test_loader
def get_dataset(path, batch_size,LM='esm2_t33_650M', ncpu=16):
    test_df = pd.read_csv(f'{path}epitope_test.tsv', sep='\t')
    train_df = pd.read_csv(f'{path}epitope_train.tsv', sep='\t')
    val_df = pd.read_csv(f'{path}epitope_val.tsv', sep='\t')
    if LM:
        test_loader = DataLoader(EpitopeDataset(test_df,f'{path}{LM}_embedding/'), batch_size=batch_size,shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        train_loader = DataLoader(EpitopeDataset(train_df,f'{path}{LM}_embedding/'), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        val_loader = DataLoader(EpitopeDataset(val_df,f'{path}{LM}_embedding/'), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
    else:
        test_loader = DataLoader(EpitopeDataset(test_df,"none"), batch_size=batch_size,shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        train_loader = DataLoader(EpitopeDataset(train_df,"none"), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        val_loader = DataLoader(EpitopeDataset(val_df,"none"), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
    return train_loader, val_loader, test_loader

def make_distance_maps(pdbfile, ca_trRosetta=False,cb_trRosetta=False,chain=None, sequence=None):
    """
    Generate (diagonalized) C_alpha and C_beta distance matrix from a pdbfile
    """
    # decide parser
    if '.gz' in pdbfile:
        pdb_handle = gzip.open(pdbfile, 'rt')
    else:
        pdb_handle = open(pdbfile,'r')
    structure_container = build_structure_container_for_pdb(pdb_handle.readlines(), chain).with_seqres(sequence)
    # structure_container.chains = {chain: structure_container.chains[chain]}
    mapper = DistanceMapBuilder(atom="CA", glycine_hack=-1)  # start with CA distances
    ca = mapper.generate_map_for_pdb(structure_container,trRosetta=ca_trRosetta)
    cb = mapper.set_atom("CB").generate_map_for_pdb(structure_container,trRosetta=cb_trRosetta)
    pdb_handle.close()
    return ca.chains, cb.chains

def retrieve_pdb(pdbfile,ca_trRosetta=False,cb_trRosetta=False, chain=None, chain_seqres=None):

    ca, cb = make_distance_maps(pdbfile,ca_trRosetta=ca_trRosetta,cb_trRosetta=cb_trRosetta,chain=chain,sequence=chain_seqres)
    assert len(ca.keys()) == 1
    chain=list(ca.keys())[0]
    angle = None
    if ca_trRosetta:
        angle = ca[chain]['trRosetta-feat']
    elif cb_trRosetta:
        angle = cb[chain]['trRosetta-feat']
    return ca[chain]['contact-map'], cb[chain]['contact-map'],angle

def write_structure_npz(ID,pdb_path,select_terms,seq,chain=None, out_dir=None):
    # there is bug about embedded_seq, since the embedded_seq is the seq itself
    """
    Write to *.npz file format.
    """
    
    A_ca, A_cb, angle = retrieve_pdb(pdbfile=pdb_path+ID,ca_trRosetta=True,chain=chain,chain_seqres=seq,cb_trRosetta=False)
    # language model embedding
    seqs = [seq]
    esm_model, alphabet = esm.pretrained.esm1_t6_43M_UR50S()  # esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    _, _, batch_tokens = batch_converter([('', x) for x in seqs])
    with torch.no_grad():
        results = esm_model(batch_tokens, repr_layers=[6], return_contacts=True)
    embedded_seq = results["representations"][6][0][1:]
    os.makedirs(out_dir, exist_ok=True)
    np.savez_compressed(os.path.join(out_dir,ID),
                        id=np.array(ID),
                        C_alpha=np.array(A_ca),
                        C_beta=np.array(A_cb),
                        trRosetta = np.array(angle),
                        labels=np.array(np.zeros(len(select_terms), dtype=np.int64)),
                        seqres=np.array(seq),
                        embedded_seq=embedded_seq
                        )

def get_dataset_by_id(path, IDs, batch_size, select_terms=None,pdb_path=None,seq=None,chain=None, ncpu=8):
    test_npzs = []

    for ID in IDs:
        if os.path.isfile(path + ID + '.npz'):
            test_npzs.append(path + ID + '.npz')
        else:
            print(f'file {ID} is not exist, and writing file to {path}')
            # if the npz file is not available, create one
            write_structure_npz(ID=ID,pdb_path=pdb_path,select_terms=select_terms,seq=seq,chain=chain, out_dir=path)
            test_npzs.append(path + ID + '.npz')

    test_loader = DataLoader(DeepfriDataset(test_npzs), batch_size=batch_size, num_workers=ncpu, pin_memory=True)

    return test_loader


def saliency_map(input_grads):
    # print('saliency_map')
    node_saliency_map = []
    for n in range(input_grads.shape[0]): # nth node
        node_grads = input_grads[n,:]
        node_saliency = torch.norm(F.relu(node_grads)).item()
        #node_saliency =torch.mean(node_grads)
        node_saliency_map.append(node_saliency)
    #node_saliency_map = MinMaxScaler(feature_range=(0, 1)).fit_transform(np.array(node_saliency_map).reshape(-1, 1)).reshape(-1, )
    return node_saliency_map

def gradcam(activation,gradient):
    node_heatmap = []
    alphas = torch.mean(gradient, axis=0)
    alphas[alphas < 0] = 0
    norm_alphas = (alphas - alphas.min()) / (alphas.max() - alphas.min())
    # print(norm_alphas.tolist())
    for n in range(activation.shape[0]):
        node_heat = F.relu(norm_alphas @ activation[n]).item()
        node_heatmap.append(node_heat)
    node_heatmap =MinMaxScaler(feature_range=(0,1)).fit_transform(np.array(node_heatmap).reshape(-1, 1)).reshape(-1, )
    df = pd.DataFrame(node_heatmap)
    return df

def write_explain_sh_script(dataset,function,label,index):
    df = pd.read_csv(dataset,sep='\t')
    zinc_df = df[df.Function.str.contains(function)]
    with open(f'slurm_explain_{label}.sh','w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -c 4                               # Request cores\n')
        f.write('#SBATCH -N 1                              # Request one node (if you request more than one core with -n, also using\n')
        f.write('                                           # -N 1 means all cores will be on the same node)\n')
        f.write('#SBATCH -p gpu\n')
        f.write('#SBATCH --gres=gpu:1\n')
        f.write('#SBATCH -t 0-02:00                         # Runtime in D-HH:MM format\n')
        f.write('#SBATCH --qos=short                       # QoS policy\n')
        f.write('#SBATCH --mem=256G                          # Memory total in MB (for all cores)\n')
        f.write('#SBATCH -o slurm_out/hostname_%j.out                 # File to which STDOUT will be written, including job ID\n')
        f.write('#SBATCH -e slurm_out/hostname_%j.err                 # File to which STDERR will be written, including job ID\n')
        f.write('#SBATCH --mail-type=END,FAIL                  # Type of email notification- BEGIN,END,FAIL,ALL\n')
        f.write('#SBATCH --mail-user=wangy549@gene.com        # Email to which notifications will be sent\n')
        f.write('hostname\n')
        f.write('ml list\n')
        f.write('ml purge\n')
        f.write('ml Anaconda3/5.0.1\n')
        f.write('ml htop\n')
        f.write('source activate GNN\n')
        f.write(f'python explain2.py -n onehot_omega -f onehot omega -lm onehot -ckn epoch=25-step=43264.ckpt -ci {index} -o /gstore/home/wangy549/deepfri/result/{label}/ -i ' + ' '.join(zinc_df.Uniprot_ID.tolist()))

