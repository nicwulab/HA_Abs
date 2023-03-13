# load dataset 

import pandas as pd
import os 
from data_preprocessing import *
from datetime import datetime

from tensorflow.keras import backend as K
from model_training import *


work_dir = '/home/yhyeung2/CoV_Encoder_HA/src/'
data_dir = work_dir + 'data/experiment2/'
cdr_char = 'XEDRKHQNSTPGCAVILMFYW-'
test_size = 0.1



# experiment 1

seq = ['CDRH1_AA', 'CDRH2_AA', 'CDRH3_AA', 'CDRL1_AA', 'CDRL2_AA', 'CDRL3_AA']
log_dir = work_dir + 'log/experiment1/'
inputx_file = 'HA_X.csv'
inputy_file = 'HA_y.csv'

# # experiment 2

# seq = ['VH_AA','VL_AA']
# log_dir = work_dir + 'log/experiment2/'
# inputx_file = 'HA-ANARCI_X.csv'
# inputy_file = 'HA-ANARCI_y.csv'



def train(X_df, y_df, random_seed, ratio):

    pad_len = []
    for c in X_df[seq]:
        pad_len.append(X_df[c].str.len().max())
    codes_dict = {i: c for i, c in enumerate(cdr_char)}
    [train_set, val_set, test_set], _ = encode(X_df, y_df, cdr_char, test_size, pad_len, le=False, col_names={'id':'Id','sequence':seq})
    [train_set_tf, val_set_tf, test_set_tf] = [to_tf_dataset(x, y, convert_dict=codes_dict) for (x,y) in [train_set, val_set, test_set]] ## debug
    (test_x_tf, y_test_tf) = test_set_tf

    (train_x, y_train), (val_x, y_val), (test_x, y_test) = [train_set, val_set, test_set]
    input_length = train_x.shape[1]



    # transformer model

    CDR_model = CDR_model_single(max_length=input_length)
    CDR_model, CDR_history = train_dl(CDR_model, (train_x, y_train), (val_x, y_val))

    model_score = CDR_model.evaluate(test_set[0], test_set[1], verbose=0)

    msg = {'time': datetime.now().strftime('%m/%d/%Y %H:%M:%S'),
           'total_params': int(CDR_model.count_params()),
           'ratio': float(ratio),
            'loss': float(model_score[0]),
            'accuracy': float(model_score[1]),
            'precision': float(model_score[2]),
            'recall': float(model_score[3]),
            'auc': float(model_score[4]),
            'prc': float(model_score[5]),
            'seed': int(random_seed),
           }

    with open(log_dir+f'single.log', 'a') as logger:
        logger.write(str(msg)+'\n')

    # dense 

    K.clear_session()
    CDR_model = CDR_model_dense(max_length=input_length)
    CDR_model, CDR_history = train_dl(CDR_model, (train_x, y_train), (val_x, y_val))

    model_score = CDR_model.evaluate(test_set[0], test_set[1], verbose=0)

    msg = {'time': datetime.now().strftime('%m/%d/%Y %H:%M:%S'),
           'total_params': int(CDR_model.count_params()),
           'ratio': float(ratio),
            'loss': float(model_score[0]),
            'accuracy': float(model_score[1]),
            'precision': float(model_score[2]),
            'recall': float(model_score[3]),
            'auc': float(model_score[4]),
            'prc': float(model_score[5]),
            'seed': int(random_seed),
           }

    with open(log_dir+f'dense.log', 'a') as logger:
        logger.write(str(msg)+'\n')

    # random forest

    K.clear_session()
    CDR_model = tfdf.keras.RandomForestModel(task=tfdf.keras.Task.CLASSIFICATION)
    CDR_model, CDR_history = train_tree(CDR_model, train_set, val_set)

    model_score = CDR_model.evaluate(test_x, y_test, verbose=0)

    msg = {'time': datetime.now().strftime('%m/%d/%Y %H:%M:%S'),
           'total_params': int(CDR_model.count_params()),
           'ratio': float(ratio),
            'loss': float(model_score[0]),
            'accuracy': float(model_score[5]),
            'precision': float(model_score[6]),
            'recall': float(model_score[7]),
            'auc': float(model_score[8]),
            'prc': float(model_score[9]),
            'seed': int(random_seed),
           }

    with open(log_dir+f'tree.log', 'a') as logger:
        logger.write(str(msg)+'\n')

def main():
    num_repeat = 6

    X_df = pd.read_csv(f'{data_dir}{inputx_file}')
    y_df = pd.read_csv(f'{data_dir}{inputy_file}')['binding']

    for m in ['single', 'dense', 'tree']:
        logger_path = log_dir+f'{m}.log'
        if not os.path.exists(os.path.dirname(logger_path)):
            os.makedirs(os.path.dirname(logger_path))
        with open(logger_path, 'a') as logger:
            logger.write('\n'*3+'*'*10 + 'START:' + datetime.now().strftime('%m/%d/%Y %H:%M:%S')+'*'*10+'\n')

    for random_seed in range(num_repeat):
        for ratio in range(1, 11):
            X_df_, y_df_ = subsample_rbd(X_df.copy(), y_df.copy(), random_seed, float(ratio)/10)
            train(X_df_, y_df_, random_seed, float(ratio)/10)

if __name__ == '__main__':
    main()