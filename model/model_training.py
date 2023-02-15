import pandas as pd
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.metrics import TruePositives, FalsePositives, TrueNegatives, FalseNegatives, BinaryAccuracy, Precision, Recall, AUC
from tensorflow.keras.metrics import CategoricalAccuracy, SparseCategoricalAccuracy
from tensorflow.keras.losses import CategoricalCrossentropy, SparseCategoricalCrossentropy
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K

from model import *

# train

def train_dl(CDR_model, train, val):
    METRICS = [
        CategoricalAccuracy(name='accuracy'),
        Precision(name='precision'),
        Recall(name='recall'),
        AUC(name='auc'),
        AUC(name='prc', curve='PR')  # precision-recall curve
    ]

    (train_x, y_train), (val_x, y_val) = train, val
    
    CDR_model.compile(
        loss=CategoricalCrossentropy(),
        optimizer=Adam(learning_rate=1e-4),
        metrics=METRICS,
    )

    callbacks = [EarlyStopping(patience=8, restore_best_weights=True)]

    CDR_history=CDR_model.fit(
        train_x, y_train,
        epochs=128, batch_size=4,
        validation_data=(val_x, y_val),
        callbacks=callbacks,
        verbose=0
    )
    return CDR_model, CDR_history

def recall_m(y_true, y_pred):
    y_true = K.ones_like(y_true)
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    all_positives = K.sum(K.round(K.clip(y_true, 0, 1)))

    recall = true_positives / (all_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    y_true = K.ones_like(y_true)
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))

    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1_score(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

def train_tree_multi(CDR_model, train, val):

    (train_x, y_train), (val_x, y_val) = train, val

    CDR_model.compile(
        metrics=[
            CategoricalAccuracy(name='accuracy'),
            f1_score, 
            precision_m, 
            recall_m,
            AUC(name='auc'),
            AUC(name='prc', curve='PR'),  # precision-recall curve
            ])

    callbacks = [EarlyStopping(patience=10, restore_best_weights=True)]

    CDR_history=CDR_model.fit(
        train_x, y_train,
        validation_data=(val_x, y_val),
        callbacks=callbacks,
        verbose=0
    )

    return CDR_model, CDR_history

def train_tree(CDR_model, train_set_tf, val_set_tf):

    (train_x, y_train), (val_x, y_val) = train_set_tf, val_set_tf

    METRICS = [
        TruePositives(name='tp'),
        FalsePositives(name='fp'),
        TrueNegatives(name='tn'),
        FalseNegatives(name='fn'),
        BinaryAccuracy(name='accuracy'),
        Precision(name='precision'),
        Recall(name='recall'),
        AUC(name='auc'),
        AUC(name='prc', curve='PR'),  # precision-recall curve
	]

    CDR_model.compile(metrics=METRICS)

    callbacks = [EarlyStopping(patience=10, restore_best_weights=True)]

    CDR_history=CDR_model.fit(
        train_x, y_train,
        validation_data=(val_x, y_val),
        callbacks=callbacks,
        verbose=0
    )

    return CDR_model, CDR_history



# test

def test_model(model, data):
    test_x, y_test, X_test_id = data
    thre = 0.5
    pred = model.predict(test_x).squeeze()
    pred_bool = pred>thre
    label = y_test#>thre
    acc = pred_bool == label
    
    eval_df = pd.DataFrame(
            {'Id':X_test_id[~acc], 
             'Label': label[~acc],
             'Prediction': pred_bool[~acc],
             'Probability':pred[~acc]
            }
        ).sort_values('Label')
#     eval_df.to_csv(f'{data_dir}/error_eval.csv')

    display(pd.Series(acc).value_counts())

    return eval_df, pred, (label, pred_bool)

def test_model_multi(model, data, le):
    test_x, y_test, X_test_id = data
    pred = model.predict(x=test_x).squeeze()
    pred_bool = pred.argmax(1)
    label = y_test.argmax(1) if len(y_test.shape)>1 else y_test
    acc = pred_bool == label
    
    eval_df = pd.DataFrame(
            {'Id':X_test_id[~acc], 
             'Label': le.inverse_transform(label[~acc]),
             'Prediction': le.inverse_transform(pred_bool[~acc]),
             'Probability Stem':pred[~acc][:,2],
             'Probability Head':pred[~acc][:,0],
             'Probability Others':pred[~acc][:,1],
            }
        ).sort_values('Prediction').sort_values('Label')

    pred_df = pd.DataFrame(
        {'Id':X_test_id, 
         'Label': le.inverse_transform(label),
         'Prediction': le.inverse_transform(pred_bool),
         'Probability Stem':pred[:,2],
         'Probability Head':pred[:,0],
         'Probability Others':pred[:,1],
        }
    ).sort_values('Prediction').sort_values('Label')
    
    display(pd.Series(acc).value_counts())
    
    return eval_df, pred, (label, pred_bool)

from sklearn.metrics import confusion_matrix
import seaborn as sns
from matplotlib import pyplot as plt

def plot_cm(labels, predictions, p=0.5, relative=False, annot=True):
    con_mat = confusion_matrix(labels, predictions > p)
    if relative:
        con_mat = np.around(con_mat.astype('float') / con_mat.sum(axis=1)[:, np.newaxis], decimals=2)
    
    plt.figure(figsize=(5,5))
    sns.heatmap(con_mat, annot=annot, fmt="d")
    plt.title('Confusion matrix p={:.2f}'.format(p))
    plt.ylabel('Actual label')
    plt.xlabel('Predicted label')
    plt.show()

def plot_cm_multi(labels, y_pred_cate, le, p=0.5, num_classes=3, relative=False, annot=True):
    con_mat = tf.math.confusion_matrix(labels=labels, predictions=y_pred_cate,num_classes=num_classes).numpy()
    
    if relative:
        con_mat = np.around(con_mat.astype('float') / con_mat.sum(axis=1)[:, np.newaxis], decimals=2)

    con_mat_df = pd.DataFrame(con_mat,
                         index = le.classes_, 
                         columns = le.classes_)
    figure = plt.figure(figsize=(5, 5))
    sns.heatmap(con_mat_df, annot=annot)#,cmap=plt.cm.Blues)
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.show()
