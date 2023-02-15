import pandas as pd

# read raw files and filter

HEADER = ['Name', 'VH_AA','VL_AA', 'CDRH1_AA','CDRH2_AA','CDRH3_AA','CDRL1_AA','CDRL2_AA','CDRL3_AA']
HEADER_FIILTER = ['CDRH1_AA','CDRH2_AA','CDRH3_AA','CDRL1_AA','CDRL2_AA','CDRL3_AA'] # filter w.r.t. CDR sequences

def read_df(original, sort_headers=None):

    '''
    headers: target columns 
    sort_headers: columns on which to sort, remove duplicates and incomplete rows
    '''
    
    original_complete = original[original[sort_headers].notnull().all(1)] # remove incomplete sequence w.r.t. sort_headers

    drop_dup = original_complete.drop_duplicates(subset=sort_headers)

    return drop_dup

def merge_columns(df, cdr_headers):

    return df.fillna('')[cdr_headers].agg(''.join, axis=1)

def concat_tables(tables, labels=None):

    if not labels:
        labels = [False, True]
        
    X_df = pd.concat(tables)
    X_df = X_df[HEADER]

    X_df['binding'] = pd.concat([pd.Series(labels[i], index=t.index) for i,t in enumerate(tables)])
    
    X_df = X_df.reset_index(drop=True)
    
    return X_df
    
def decouple_headstemothers(neg_excel, pos_excel_head, pos_excel_stem):

    negatives = read_df(neg_excel, HEADER_FIILTER)
    positives_head = read_df(pos_excel_head, HEADER_FIILTER) # all +ve cases in S_HA
    positives_stem = read_df(pos_excel_stem, HEADER_FIILTER)
    
    # process duplicates

    negative_cols = merge_columns(negatives, HEADER_FIILTER)
    positive_cols_h = merge_columns(positives_head, HEADER_FIILTER)
    positive_cols_s = merge_columns(positives_stem, HEADER_FIILTER)
    
    # ensure all HA contained in all_paired

    negative_mask = negative_cols.isin(positive_cols_h) & negative_cols.isin(positive_cols_s)
    negatives = negatives.loc[~negative_mask] # all other -ve cases; removed rows duplicated in HA / +ve

    return negatives, positives_head, positives_stem

def process(All, HAhead, HAstem):

    neg, pos_head, pos_stem = decouple_headstemothers(All, HAhead, HAstem)

    neg_mask = neg['Name'].isin(pos_head['Name']) | neg['Name'].isin(pos_stem['Name'])
    neg = neg.loc[~neg_mask]

    X_df = concat_tables([pos_stem, pos_head, neg], ['stem', 'head', 'others'])

    X_df = read_df(X_df, ['VH_AA','VL_AA'])
    y_df = X_df[['Name','binding']]
    X_df = X_df.drop('binding',axis=1)
    
    return X_df, y_df

# ANARCI sequence alignment

def process_duplicate(df, index_col, suffix_col):

    df[index_col] = df[index_col].where(~df[index_col].duplicated(keep=False),df[index_col].str.cat(df[suffix_col].astype(str), sep='_'))

def anarci_process(x_h, x_l, y_original=None):
    
    # x_h_dup_index = x_h['Id'].duplicated(keep=False)
    # x_l_dup_index = x_l['Id'].duplicated(keep=False)
    
    process_duplicate(x_h, 'Id', 'seqend_index')
    process_duplicate(x_l, 'Id', 'seqend_index')

    x_h_ = x_h.iloc[:, 13:].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    x_l_ = x_l.iloc[:, 13:].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

    x_h_ = pd.DataFrame({'VH_AA': x_h_, 'Id':x_h['Id']})
    x_l_ = pd.DataFrame({'VL_AA': x_l_, 'Id':x_l['Id']})

    x = pd.DataFrame(x_l_).merge(x_h_, how='inner', on='Id')
    x = x[['Id', 'VH_AA', 'VL_AA']]

    xh_in_xl = x_h['Id'].isin(x_l['Id'])
    xl_in_xh = x_l['Id'].isin(x_h['Id'])
    
    join_rows = x_h.loc[xh_in_xl,'Id']
    merge_rows = x['Id']
    equal =  join_rows.reset_index(drop=True).equals(merge_rows.reset_index(drop=True))
    
    if equal:
        x = x.reset_index(drop=True)
    else:
        # display('join', x_h.loc[xh_in_xl,'Id'])
        # display('merge', x['Id'])
        raise ValueError('Duplicated sequences found')
    
    return x

def decouple_haothers(neg_excel, pos_excel):
    
    negatives = read_df(neg_excel, HEADER_FIILTER)
    positives = read_df(pos_excel, HEADER_FIILTER) # all +ve cases in S_HA
    
    # process duplicates

    negative_cols = merge_columns(negatives, HEADER_FIILTER)
    positive_cols = merge_columns(positives, HEADER_FIILTER)
    
    # ensure all HA contained in all_paired

    negative_mask = negative_cols.isin(positive_cols)
    negatives = negatives.loc[~negative_mask] # all other -ve cases; removed rows duplicated in HA / +ve

    X_df = pd.concat([negatives, positives])
    X_df = X_df[HEADER]

    binding = pd.concat([pd.Series(False, index=negatives.index), pd.Series(True, index=positives.index)])
    y_df = pd.DataFrame({'Name': X_df['Name'], 'binding': binding})
    
    X_df = X_df.reset_index(drop=True)
    y_df = y_df.reset_index(drop=True)
    
    return X_df, y_df

# encoding

from sklearn.utils import shuffle
from keras_preprocessing.sequence import pad_sequences
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import tensorflow_decision_forests as tfdf

def subsample_stemhead(X_df, y_df, random_seed=0, neg_ratio=1):
    
    all_paired_mask = y_df == 'others'

    cleaned_df_HA = X_df.loc[~all_paired_mask]
    cleaned_df_HA_y = y_df.loc[~all_paired_mask]
    
    cleaned_df = X_df.loc[all_paired_mask]
    
#     assuming classes = {'stem', 'head', 'others'}, num_sample of -ve cases being average of 'stem' and 'head'
    num_sample = int((~all_paired_mask).sum()/2) * neg_ratio

    # false cases
    cleaned_df = shuffle(cleaned_df, random_state=random_seed, n_samples=num_sample)

    X_df = pd.concat([cleaned_df, cleaned_df_HA])
    y_df = pd.concat([pd.Series('others', index=range(num_sample)), cleaned_df_HA_y])
    
    X_df = X_df.reset_index(drop=True)
    y_df = y_df.reset_index(drop=True)
    
    return X_df, y_df

def subsample_rbd(X_df, y_df, random_seed=0, ratio=1.0):
    
    rbd_mask = y_df == False

    cleaned_df_HA = X_df.loc[~rbd_mask]
    cleaned_df = X_df.loc[rbd_mask]
    
#     assuming classes = {'stem', 'head', 'others'}, num_sample of -ve cases being average of 'stem' and 'head'
    num_sample = int((~rbd_mask).sum()*ratio)
    cleaned_df_HA = shuffle(cleaned_df_HA, random_state=random_seed, n_samples=num_sample)

    # false cases
    cleaned_df = shuffle(cleaned_df, random_state=random_seed, n_samples=num_sample)

    X_df = pd.concat([cleaned_df, cleaned_df_HA])
    y_df = pd.concat([pd.Series(False, index=range(num_sample)), pd.Series(True, index=range(num_sample))])
    
    X_df = X_df.reset_index(drop=True)
    y_df = y_df.reset_index(drop=True)

    return X_df, y_df


def integer_encoding(df, cdr_char):
    '''
    cdr_char = 'XEDRKHQNSTPGCAVILMFYW-'
    'X': null character for padding masking
    '-': alignment space character for position information
    '''

    if df.isnull().values.any():
        raise RuntimeError('DataFrame contains empty/ nan cells...')
    codes = list(cdr_char)

    # unknown = 0
    return df.applymap(lambda x: list(map(lambda z: codes.index(z) if z in codes else 0, x)))


def padding(df, padding_maxlen):
    col_list = []
    for i in range(df.shape[1]):
        col_list.append(pad_sequences(df[:,i], maxlen=padding_maxlen[i], padding='post', truncating='post'))
    return col_list

def encode_six_CDR(data, cdr_char, padding_maxlen):
    """
    - Encodes code sequence to integer values.
    - add post-padding
    """

    data_integer = integer_encoding(data, cdr_char)
    data_pad = padding(data_integer.to_numpy(), padding_maxlen)

    return data_pad


def encode(X_df, y_df, cdr_char, test_size, pad=30, test_only=False, le=True, col_names={'id':'Id','sequence':['VH_AA','VL_AA']}):
    
    if not test_only:
        X_train, X_test, y_train, y_test = train_test_split(X_df, y_df, test_size=test_size)
        X_train, X_eval, y_train, y_eval = train_test_split(X_train, y_train, test_size=test_size)
    else:
        (X_train, X_eval, X_test), (y_train, y_eval, y_test) = X_df, y_df

    col_id = col_names['id']
    col_seq = col_names['sequence']
    X_train_id, X_eval_id, X_test_id = X_train[col_id], X_eval[col_id], X_test[col_id]
    X_train, X_eval, X_test = X_train[col_seq], X_eval[col_seq], X_test[col_seq]

    train_x = encode_six_CDR(X_train, cdr_char, pad)
    val_x = encode_six_CDR(X_eval, cdr_char, pad)
    test_x = encode_six_CDR(X_test, cdr_char, pad)

    train_x = merge_cdr(train_x); val_x = merge_cdr(val_x); test_x = merge_cdr(test_x)
    
    # LabelEncoder
    if le:
        le = LabelEncoder()

        y_train = le.fit_transform(y_train)
        y_eval = le.transform(y_eval)
        y_test = le.transform(y_test)
        
        return [(train_x, y_train), (val_x, y_eval), (test_x, y_test)], [X_train_id, X_eval_id, X_test_id], le
    
    return [(train_x, y_train), (val_x, y_eval), (test_x, y_test)], [X_train_id, X_eval_id, X_test_id]

def merge_cdr(cdr_list):
    df = pd.DataFrame(cdr_list[0])
    for cdr in cdr_list[1:]:
        df = pd.concat([df, pd.DataFrame(cdr)], axis=1)
    
    return df.to_numpy()
    
def to_tf_dataset(x, y, y_label='binding',convert_dict=None):
    x = pd.DataFrame(x)
    if convert_dict:
        x = x.replace(convert_dict)
    y = pd.Series(y).reset_index(drop=True)
    return x, y
