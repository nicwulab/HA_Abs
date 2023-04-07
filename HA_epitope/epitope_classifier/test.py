import torch
import pytorch_lightning as pl
import numpy as np
import pandas as pd
import os
import argparse
from model import Epitope_Clsfr
from utils import get_dataset
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, multilabel_confusion_matrix
# Usage
'''
python test.py
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batch_size", default=32, type=int)
    parser.add_argument("-lr", "--learning_rate", default=1e-4, type=float)
    parser.add_argument("-c", "--classes", default=7, type=int)
    parser.add_argument("-lm", "--language_model", default='esm2_t33_650M_UR50D', type=str)
    parser.add_argument("-l", "--layers", default=3, type=int)
    parser.add_argument("-hd", "--hidden_dim", default=1280, type=int)
    parser.add_argument("-dp", "--dataset_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/result/', type=str)
    parser.add_argument("-ckp", "--checkpoint_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/checkpoint/', type=str) 
    parser.add_argument("-ckn","--checkpoint_name", default='epoch=205-step=73954.ckpt', type=str)
    parser.add_argument("-n", "--name", default='esm_baseline_new', type=str)
    parser.add_argument("-o", "--output_path", default='/home/yiquan2/ESM_Ab/Ab_epitope/result/', type=str)
    args = parser.parse_args()


    train_loader, val_loader, test_loader = get_dataset(args.dataset_path, batch_size=args.batch_size)
    trainer = pl.Trainer()
    # Check whether pretrained model exists. If yes, load it and skip training
    pretrained_filename = os.path.join(args.checkpoint_path+args.name+'/', args.checkpoint_name)
    if os.path.isfile(pretrained_filename):
        print("Found pretrained model, loading...")
        model = Epitope_Clsfr.load_from_checkpoint(pretrained_filename,classes=args.classes,hidden_dim=args.hidden_dim,layers=args.layers,class_weights=None,lm_model_name = args.language_model)
        model.eval()
        # test on test set 
        test_result = trainer.test(model,test_loader)
        # test on val set 
        val_result = trainer.test(model,val_loader)

    # export prediction table and confusion matrix
        # create empty lists to store the true and predicted labels
        true_labels_ls = []
        predicted_labels_ls = []
        predicted_probabilities = []

        classes = ["HA:Head", "HA:Stem","HIV", "S:NTD", "S:RBD", "S:S2", "Others"]
        # loop over the test data and predict the labels
        with torch.no_grad():
            for batch in test_loader:
                # get the inputs and labels
                inputs, labels = batch
                outputs = model(inputs)
                # get model prediction probabilities
                probs = torch.softmax(outputs,dim=1)
                # get model predicted classes index
                predicted = torch.argmax(outputs, dim=1)

                predicted_probabilities.append(probs)
                true_labels_ls.append(labels)
                predicted_labels_ls.append(predicted)

        labels_all = np.concatenate(true_labels_ls)
        predicted_all = np.concatenate(predicted_labels_ls)
        probabilities_all = np.concatenate(predicted_probabilities)
        # read df
        df = pd.read_csv(f'{args.output_path}epitope_test.tsv',sep='\t')
        df['predicted_class'] = predicted_all
        df['real_label'] = labels_all
        df['predicted_probability'] = probabilities_all[np.arange(len(predicted_all)), predicted_all]
        df.to_csv(f'{args.output_path}{args.name}_epitope_test_prediction.tsv',sep='\t')

        # create the confusion matrix
        cm = confusion_matrix(labels_all, predicted_all)
        # plot the confusion matrix
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        # classes = [str(i) for i in range(cm.shape[0])]
        plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)

        plt.title("Normalized confusion matrix")
        plt.colorbar()
        tick_marks = np.arange(len(classes))
        plt.xticks(tick_marks, classes, rotation=45)
        plt.yticks(tick_marks, classes)
        plt.tight_layout()
        plt.ylabel('True label')
        plt.xlabel('Predicted label')
        plt.savefig(f'{args.output_path}{args.name}_confusion_matrix.png', dpi=300)





