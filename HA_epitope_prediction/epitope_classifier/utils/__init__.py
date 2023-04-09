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


def load_esm_embedding(path):

    esm_embedding = torch.load(path)['representations'][33]
    return esm_embedding


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
        self.label_map={"HA:Head":0, "HA:Stem":1, "Others":2}
        self.load_esm_embedding = partial(load_esm_embedding)
    def __len__(self):
         return len(self.df)
    def __getitem__(self, index):
        row = self.df.iloc[index]
        name = row['Name']
        seq = row['VH_AA']
        target_data = row['Antigen_epitopes']
        # Convert your input and target data into PyTorch tensors
        if self.embed_path != "none":
            esm_embedded_tensor = self.load_esm_embedding(f'{self.embed_path}{name}.pt')
            x=zero_padding(esm_embedded_tensor,max_h=150,max_w=1280)
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
    fas_f = f'/home/yiquan2/ESM_Ab/HA_Abs/HA_epitope_prediction/result/{filename}.fasta'
    if not os.path.exists(fas_f):
        print('start generating fasta file')
        df2fasta(df,fas_f)
        # remove duplicate sequence from fasta
        subprocess.run(['python', '/home/yiquan2/ESM_Ab/HA_Abs/HA_epitope_prediction/rm_fas_repeats.py', fas_f, fas_f])
    # run esm embedding
    esm_output_path = f'/home/yiquan2/ESM_Ab/HA_Abs/HA_epitope_prediction/result/{filename}/'
    os.makedirs(esm_output_path, exist_ok=True)
    if len(os.listdir(esm_output_path)) < 1:
        print('start embedding seqeuence...')
        subprocess.run(['python', '/home/yiquan2/ESM_Ab/HA_Abs/HA_epitope_prediction/esm_extractor.py',  'esm2_t33_650M_UR50D',fas_f, f'{esm_output_path}', '--repr_layers', '33', '--include', 'per_tok','contacts'])
    test_loader = DataLoader(EpitopeDataset(df,esm_output_path,prediction=True), batch_size=batch_size,shuffle=False,
                                         num_workers=ncpu, pin_memory=True)
    return test_loader
def get_dataset(path, batch_size,LM=True, ncpu=16):
    test_df = pd.read_csv(f'{path}epitope_test.tsv', sep='\t')
    train_df = pd.read_csv(f'{path}epitope_train.tsv', sep='\t')
    val_df = pd.read_csv(f'{path}epitope_val.tsv', sep='\t')
    if LM:
        test_loader = DataLoader(EpitopeDataset(test_df,f'{path}esm2_t33_650M_embedding/'), batch_size=batch_size,shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        train_loader = DataLoader(EpitopeDataset(train_df,f'{path}esm2_t33_650M_embedding/'), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        val_loader = DataLoader(EpitopeDataset(val_df,f'{path}esm2_t33_650M_embedding/'), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
    else:
        test_loader = DataLoader(EpitopeDataset(test_df,"none"), batch_size=batch_size,shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        train_loader = DataLoader(EpitopeDataset(train_df,"none"), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
        val_loader = DataLoader(EpitopeDataset(val_df,"none"), batch_size=batch_size, shuffle=False,
                                num_workers=ncpu, pin_memory=True)
    return train_loader, val_loader, test_loader

