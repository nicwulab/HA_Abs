import torch
import pytorch_lightning as pl
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from pytorch_lightning.callbacks.model_checkpoint import ModelCheckpoint
from pytorch_lightning.loggers import WandbLogger
from sklearn.utils.class_weight import compute_class_weight
import numpy as np
import os
import argparse
from utils import get_dataset
from model import Epitope_Clsfr

# Usage
'''
python train.py
'''

def train_Epitope_Clsfr(model_name,batch_size,classes,lm_model_name,lr,dataset_path,CHECKPOINT_PATH,ckn,hidden_dim=1280,layers=3,class_weights=None,logger=None):
    
    # dataloader
    if  lm_model_name == 'onehot':
        train_loader, val_loader, test_loader = get_dataset(dataset_path, batch_size=batch_size,LM=False)
    elif  lm_model_name == 'mBLM':
        train_loader, test_loader, val_loader = get_dataset(dataset_path, batch_size=batch_size,LM='mBLM')
    elif  lm_model_name == 'esm2_t33_650M_UR50D':
        train_loader, val_loader, test_loader = get_dataset(dataset_path, batch_size=batch_size,LM='esm2_t33_650M')
    else:
        train_loader, val_loader, test_loader = get_dataset(dataset_path, batch_size=batch_size)
    pl.seed_everything(42)
    # Create a PyTorch Lightning trainer with the generation callback
    root_dir = os.path.join(CHECKPOINT_PATH, model_name)
    os.makedirs(root_dir, exist_ok=True)
    

    trainer = pl.Trainer(default_root_dir=root_dir,
                         logger=logger,
                         max_epochs=25,
                         strategy='ddp_find_unused_parameters_false',
                         callbacks=[ModelCheckpoint(dirpath=root_dir, filename='{epoch}_'+ckn, every_n_epochs = 1, save_top_k = -1, save_weights_only=True),
                                    EarlyStopping(monitor="val_loss", mode="min", patience=5)])

    model = Epitope_Clsfr(classes=classes,class_weights=class_weights,lr=lr,lm_model_name = lm_model_name,hidden_dim=hidden_dim,layers=layers)
    trainer.fit(model, train_loader, val_loader)
    best_model_path = trainer.checkpoint_callback.best_model_path
    print('best model path is ', best_model_path)
    model = Epitope_Clsfr.load_from_checkpoint(trainer.checkpoint_callback.best_model_path,classes=classes,class_weights=class_weights,lm_model_name = lm_model_name,hidden_dim=hidden_dim,layers=layers)
    # Test best model on validation and test set
    # val_result = trainer.test(model, val_loader, verbose=False)
    test_result = trainer.test(model, test_loader, verbose=False)
    # result = {"test": test_result[0]['test_F1'], "val": val_result[0]['test_F1']}
    result = ''
    return model, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batch_size", default=32, type=int)
    parser.add_argument("-lr", "--learning_rate", default=2e-5, type=float)
    parser.add_argument("-c", "--classes", default=7, type=int)
    parser.add_argument("-l", "--layers", default=1, type=int)
    parser.add_argument("-hd", "--hidden_dim", default=768, type=int)
    parser.add_argument("-lm", "--language_model", default='mBLM', type=str)
    parser.add_argument("-dp", "--dataset_path", default='result/', type=str)
    parser.add_argument("-ckp", "--checkpoint_path", default='checkpoint/', type=str)
    parser.add_argument("-n", "--name", default='mBLM_attention', type=str)
    parser.add_argument("-ckn", "--checkpoint_name", default='mBLM', type=str)
    args = parser.parse_args()

    logger = WandbLogger(name=args.name, project='Epitope_clsfr')
    model, result = train_Epitope_Clsfr(model_name=args.name,
                                   classes=args.classes,
                                   lr=args.learning_rate,
                                   batch_size=args.batch_size,
                                   logger=logger,
                                   class_weights=None,
                                   lm_model_name=args.language_model,
                                   dataset_path=args.dataset_path,
                                   hidden_dim=args.hidden_dim,
                                   layers=args.layers,
                                   CHECKPOINT_PATH=args.checkpoint_path,
                                   ckn = args.checkpoint_name
                                   )
    

