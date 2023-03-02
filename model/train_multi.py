# load dataset 

import pandas as pd
import os 
from data_preprocessing import *
from datetime import datetime

from tensorflow.keras.utils import to_categorical
from tensorflow.keras import backend as K
from model_training import *

work_dir = '/home/yhyeung2/CoV_Encoder_HA/src/'
seq = ['VH_AA','VL_AA']
cdr_char = 'XEDRKHQNSTPGCAVILMFYW-'
test_size = 0.1

data_dir = work_dir + 'data/experiment3/' # 'log/experiment4/'
log_dir = work_dir + 'log/experiment3/' # 'log/experiment4/'

def train(X_df, y_df, random_seed, ratio):

    pad_len = []
    for c in X_df[seq]:
        pad_len.append(X_df[c].str.len().max())
    codes_dict = {i: c for i, c in enumerate(cdr_char)}
    [train_set, val_set, test_set], _, le = encode(X_df, y_df, cdr_char, test_size, pad_len, le=True, col_names={'id':'Id','sequence':seq})

    (train_x, y_train), (val_x, y_val), (test_x, y_test) = [train_set, val_set, test_set]
    [(val_x, y_val), (test_x, y_test)] = [sample_equal_ratio(x,y,le) for (x, y) in [(val_x, y_val), (test_x, y_test)]]
    val_set, test_set = (val_x, y_val), (test_x, y_test)
    
    y_train, y_val, y_test = to_categorical(y_train), to_categorical(y_val), to_categorical(y_test)
    (test_x, y_test) = test_set
    input_length = train_set[0].shape[1]
    
    [train_set_tf, val_set_tf, test_set_tf] = [to_tf_dataset(x.copy(), y.copy(), convert_dict=codes_dict) for (x,y) in [train_set, val_set, test_set]]
    (test_x_tf, y_test_tf) = test_set_tf
    y_test_tf = to_categorical(y_test_tf)

    # transformer model

    CDR_model = CDR_model_single(max_length=input_length, n_classes=3)
    CDR_model, CDR_history = train_dl_multi(CDR_model, (train_x, y_train), (val_x, y_val))

    model_score = CDR_model.evaluate(test_x, y_test, verbose=0)

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
    CDR_model = CDR_model_dense(max_length=input_length, n_classes=3)
    CDR_model, CDR_history = train_dl_multi(CDR_model,(train_x, y_train), (val_x, y_val))

    model_score = CDR_model.evaluate(test_x, y_test, verbose=0)

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
    CDR_model, CDR_history = train_tree_multi(CDR_model, train_set, val_set)

    model_score = CDR_model.evaluate(test_x, y_test, verbose=0)

    msg = {'time': datetime.now().strftime('%m/%d/%Y %H:%M:%S'),
           'total_params': int(CDR_model.count_params()),
           'ratio': float(ratio),
            'loss': float(model_score[0]),
            'accuracy': float(model_score[1]),
            'f1score': float(model_score[2]),
            'precision': float(model_score[3]),
            'recall': float(model_score[4]),
            'auc': float(model_score[5]),
            'prc': float(model_score[6]),
            'seed': int(random_seed),
           }

    with open(log_dir+f'tree.log', 'a') as logger:
        logger.write(str(msg)+'\n')

def main():

    X_df = pd.read_csv(f'{data_dir}HA-ANARCI_X.csv')
    y_df = pd.read_csv(f'{data_dir}HA-ANARCI_y.csv')['binding']

    false_cnt = (y_df == 'others').sum()
    true_false_ratio = false_cnt // (y_df.shape[0] - false_cnt)

    for m in ['single', 'dense', 'tree']:
        logger_path = log_dir+f'{m}.log'
        if not os.path.exists(os.path.dirname(logger_path)):
            os.makedirs(os.path.dirname(logger_path))
        with open(logger_path, 'a') as logger:
            logger.write('\n'*3+'*'*10 + 'START:' + datetime.now().strftime('%m/%d/%Y %H:%M:%S')+'*'*10+'\n')

    for random_seed in range(6):
        for ratio in range(1, true_false_ratio+1, true_false_ratio//14):
            X_df_, y_df_ = subsample_stemhead(X_df.copy(), y_df.copy(), random_seed, ratio)
            train(X_df_, y_df_, random_seed, ratio)

if __name__ == '__main__':
    main()