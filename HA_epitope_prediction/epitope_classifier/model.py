import torch
from torch import nn
from torch.nn import functional as F
from torchmetrics import Accuracy,AveragePrecision,AUROC,PrecisionRecallCurve,F1Score
import pytorch_lightning as pl



class Epitope_Clsfr(pl.LightningModule):
    """ Class containig the LM models for predicting anitbody epitope. """
    def __init__(self, classes,class_weights=None,seq_dim=26,hidden_dim=1280,lr=0.0001, drop=0.1, lm_model_name='esm2_t33_650M_UR50D',layers=2):
        super().__init__()
        """ Initialize the model
        :param classes: {int} number of epitopes
        :param class_weights: {tensor} class weights
        :param seq_dim: {int} sequence features per residue (26 for 1-hot encoding or language model embedding dimension)
        :param hidden_dim: {int} hidden features dimension (default maximum 1500)
        :param lr: {float} learning rate for Adam optimizer
        :param drop: {float} dropout fraction for Dense layers
        :lm_model: {string} name of the pre-trained ESM language model to be loaded
        """

        self.lr = lr
        self.class_weights =class_weights #.to(DEVICE)
        self.drop = drop
        self.classes = classes
        self.lm_model_name = lm_model_name
        self.layers = layers
        # LM layer
        if lm_model_name == "esm2_t33_650M_UR50D":
            self.seq_dim = 1280
            seq_dim = 1280
        elif lm_model_name == "esm2_t36_3B_UR50D":
            self.seq_dim = 2560
            seq_dim = 2560
        elif lm_model_name == "esm2_t6_8M_UR50D":
            self.seq_dim = 320
            seq_dim = 320
        elif lm_model_name == 'onehot':
            self.seq_dim = 26
            seq_dim = 26
        self.accuracy = Accuracy(task='multiclass',num_classes=classes)
        self.AveragePrecision = AveragePrecision(task='multiclass',num_classes=classes)
        self.AUCROC = AUROC(task='multiclass',num_classes=classes)
        self.PrecisionRecallCurve = PrecisionRecallCurve(task='multiclass',num_classes=classes)
        self.F1Score=F1Score(task='multiclass',num_classes=classes)
        # define classifier layers
        self.avgpool = nn.AdaptiveAvgPool2d((1,hidden_dim))
        self.fc = nn.Linear(hidden_dim, hidden_dim)

        # final layer
        self.fc2 = nn.Linear(hidden_dim, classes)

    def forward(self,x):
            
        # Sum pooling
        x = self.avgpool(x)
        x = x.view(x.size(0), -1)
        x= F.dropout(x, p=self.drop, training=self.training)
        for n in range(self.layers):
            x = self.fc(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.drop, training=self.training)
        y = self.fc2(x)
        return y

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        return optimizer

    def training_step(self, batch, batch_idx):
        x,y = batch

        x_out = self.forward(x)
        loss = F.cross_entropy(x_out, y)
        self.accuracy(x_out, y)
        self.F1Score(x_out, y)
        self.log('train_loss', loss, batch_size=y.shape[0])
        self.log('train_accuracy', self.accuracy, batch_size=y.shape[0])
        self.log('train_F1', self.F1Score, batch_size=y.shape[0])
        
        return loss

    def validation_step(self, batch, batch_index):

        x,y = batch
        x_out = self.forward(x)
        loss = F.cross_entropy(x_out, y)
        self.accuracy(x_out, y)
        self.F1Score(x_out, y)
        self.log('val_loss', loss, batch_size=y.shape[0])
        self.log('val_accuracy', self.accuracy, batch_size=y.shape[0])
        self.log('val_F1', self.F1Score, batch_size=y.shape[0])
    def test_step(self, batch, batch_index):
        x,y = batch
        x_out = self.forward(x)
        loss = F.cross_entropy(x_out, y)
        self.accuracy(x_out, y)
        self.F1Score(x_out, y)
        self.log('test_loss', loss,batch_size=y.shape[0])
        self.log('test_accuracy', self.accuracy, batch_size=y.shape[0])
        self.log('test_F1', self.F1Score, batch_size=y.shape[0])


