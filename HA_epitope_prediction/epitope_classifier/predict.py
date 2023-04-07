import torch
import pytorch_lightning as pl
import numpy as np
import pandas as pd
import os
import argparse
from model import Epitope_Clsfr
from utils import get_dataset, get_dataset_from_df
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, multilabel_confusion_matrix
# Usage
'''
python predict.py
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batch_size", default=32, type=int) 
    parser.add_argument("-lr", "--learning_rate", default=1e-4, type=float)
    parser.add_argument("-c", "--classes", default=7, type=int)
    parser.add_argument("-lm", "--language_model", default='esm2_t33_650M_UR50D', type=str)
    parser.add_argument("-l", "--layers", default=3, type=int)
    parser.add_argument("-hd", "--hidden_dim", default=1280, type=int)
    parser.add_argument("-dp", "--dataframe_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/result/Flu_unknown.csv', type=str)
    parser.add_argument("-ckp", "--checkpoint_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/checkpoint/', type=str)
    parser.add_argument("-ckn","--checkpoint_name", default='epoch=205-step=73954.ckpt', type=str)
    parser.add_argument("-n", "--name", default='esm_baseline_new', type=str)
    parser.add_argument("-o", "--output_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/result/', type=str)
    args = parser.parse_args()

    test_loader = get_dataset_from_df(args.dataframe_path, batch_size=args.batch_size)
    filename = args.dataframe_path.split('/')[-1].split('.')[0]
    pretrained_filename = os.path.join(args.checkpoint_path+args.name+'/', args.checkpoint_name)
    if os.path.isfile(pretrained_filename):
        print("Found pretrained model, loading...")
        model = Epitope_Clsfr.load_from_checkpoint(pretrained_filename,classes=args.classes,hidden_dim=args.hidden_dim,layers=args.layers,class_weights=None,lm_model_name = args.language_model)
        model.eval()
        # test on test set
        predicted_labels_ls = []
        predicted_probabilities = []

        with torch.no_grad():
            for batch in test_loader:
                # get the inputs and labels
                inputs, _ = batch

                outputs = model(inputs)
                # get model prediction probabilities
                probs = torch.softmax(outputs,dim=1)
                # get model predicted classes index
                predicted = torch.argmax(outputs, dim=1)

                predicted_labels_ls.append(predicted)
                predicted_probabilities.append(probs)

        predicted_all = np.concatenate(predicted_labels_ls)
        probabilities_all = np.concatenate(predicted_probabilities)

        # read df
        if args.dataframe_path.split('.')[-1].lower() == 'csv':
            df = pd.read_csv(args.dataframe_path)

        elif args.dataframe_path.split('.')[-1].lower() == 'tsv':
            df = pd.read_csv(args.dataframe_path,sep='\t')
        elif args.dataframe_path.split('.')[-1].lower() == 'xlsx':
            df = pd.read_excel(args.dataframe_path)
        else:
            print(f"Error: unsupported file type for {args.dataframe_pat}")
        # add the predicted class and probability to the DataFrame

        df['predicted_class'] = predicted_all
        df['predicted_probability'] = probabilities_all[np.arange(len(predicted_all)), predicted_all]
        df.to_csv(f'{args.output_path}{filename}_prediction.tsv',sep='\t')

