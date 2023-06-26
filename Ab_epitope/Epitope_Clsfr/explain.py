import torch
import pytorch_lightning as pl
import numpy as np
import pandas as pd
import os
import argparse
from model import Epitope_Clsfr
from utils import get_dataset,get_dataset_from_df, saliency_map,gradcam
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, multilabel_confusion_matrix
# Usage
'''
python explain.py
'''
activation ={}
def get_gradient(name):
    def hook(model, grad_input, grad_output):
        activation[name] = grad_output
    return hook
def get_activation(name):
    def hook(model, input, output):
        activation[name] = output.detach()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batch_size", default=1, type=int)
    parser.add_argument("-lr", "--learning_rate", default=2e-5, type=float)
    parser.add_argument("-c", "--classes", default=7, type=int)
    parser.add_argument("-lm", "--language_model", default='mBLM', type=str)
    parser.add_argument("-l", "--layers", default=3, type=int)
    parser.add_argument("-hd", "--hidden_dim", default=768, type=int)
    parser.add_argument("-dp", "--dataset_path", default='result/', type=str)
    parser.add_argument("-dfp", "--dataframe_path", default='result/Flu_unknown.csv', type=str)
    parser.add_argument( "--provide_dataset", action='store_true')
    parser.add_argument("-ckp", "--checkpoint_path", default='checkpoint/', type=str)
    parser.add_argument("-ckn","--checkpoint_name", default='epoch=19-step=7160.ckpt', type=str)
    parser.add_argument("-n", "--name", default='mBLM_attention', type=str)
    parser.add_argument("-o", "--output_path", default='result/explain_mBLM/', type=str)
    args = parser.parse_args()

    # default dataset is the training test val set
    if args.provide_dataset:
        test_loader = get_dataset_from_df(args.dataframe_path, batch_size=args.batch_size)

    else:
        train_loader, val_loader, test_loader = get_dataset(args.dataset_path, batch_size=args.batch_size,LM='mBLM')
    # Check whether pretrained model exists. If yes, load it and skip training
    pretrained_filename = os.path.join(args.checkpoint_path+args.name+'/', args.checkpoint_name)
    if os.path.isfile(pretrained_filename):
        print("Found pretrained model, loading...")
        model = Epitope_Clsfr.load_from_checkpoint(pretrained_filename,classes=args.classes,hidden_dim=args.hidden_dim,layers=args.layers,class_weights=None,lm_model_name = args.language_model)
        model.eval()

#         def backward_hook(module, grad_input, grad_output):
#             # Record the gradient of the input and output tensors
#             print("Gradient of input:", grad_input[0])
#             print("Gradient of output:", grad_output[0].size())

#         def forward_hook(module, input, output):
#             # Record the activation of the input tensor
#             print("Activation of input:", input)
#         model.attention.register_backward_hook(backward_hook)
#         model.attention.register_forward_hook(forward_hook)
        os.makedirs(args.output_path, exist_ok=True)
    # read df
        if args.provide_dataset:
            if args.dataframe_path.split('.')[-1] == 'tsv':
                df = pd.read_csv(args.dataframe_path,sep='\t')
            elif args.dataframe_path.split('.')[-1] == 'csv':
                df = pd.read_csv(args.dataframe_path)
            elif args.dataframe_path.split('.')[-1] == 'xlsx':
                df = pd.read_excel(args.dataframe_path)
            else:
                print('ERROR: the format of dataframe is not defined')
        else:    
            df = pd.read_csv(f'{args.dataset_path}epitope_train.tsv',sep='\t')
        
        names = df.loc[:,'Name']
        # loop over the test data and predict the labels

        for idx,batch in enumerate(train_loader):
            # if idx == 2:break
            # get the inputs and labels
            inputs, labels = batch

            for p in model.parameters():
                p.requires_grad = False

            inputs.requires_grad = True

            outputs = model(inputs)

            outputs[:,labels].sum().backward()
            seq_saliency_map = saliency_map(inputs.grad[0])
            gradcam_map = gradcam(inputs[0],inputs.grad[0])
            name = names[idx].replace('/', '_')
            seq_saliency_df = pd.DataFrame(seq_saliency_map)
            seq_saliency_df.to_csv(f'{args.output_path}{name}_grad.tsv',sep='\t')
            gradcam_map.to_csv(f'{args.output_path}{name}_gradcam.tsv',sep='\t')