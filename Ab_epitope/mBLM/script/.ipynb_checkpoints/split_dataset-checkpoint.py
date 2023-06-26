import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
from sklearn.model_selection import GroupShuffleSplit,train_test_split
'''
split into train~val/test(50%, 60%, 70%, 80%, 90% from training set) by different sequence identity
'''


# turn off the SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def data_split(df,clusters,out_path):
    '''
    split data into train~val/test
    :param df:
    :param out_path:
    :return:
    '''
    print("start data splitting \n")
    # create a test set list based on the sequence identity
    testset_ls = []
    for c in clusters:
        splitter = GroupShuffleSplit(test_size=.01, n_splits=2, random_state=42)
        split = splitter.split(df, groups=df[f'cluster_ID_by_{c}'])
        train_inds, test_inds = next(split)
        train_df = df.iloc[train_inds]
        test_df = df.iloc[test_inds]
        testset_ls.append(test_df)
        # update new df to split
        df = train_df
    
    # use 70% sequence identity to split training set and val set
    splitter = GroupShuffleSplit(test_size=.05, n_splits=2, random_state=42)
    split = splitter.split(df, groups=df[f'cluster_ID_by_{clusters[2]}'])
    train_inds, val_inds = next(split)
    train_df = df.iloc[train_inds]
    val_df = df.iloc[val_inds]
    
    #upsampling for the training set
    # upsample the training set
    label_counts = train_df['v_call_light'].str.split('*').str[0].value_counts()
    minority_labels = label_counts[label_counts < 5000].index.tolist()

    df_resampled = pd.DataFrame(columns=train_df.columns)

    for label in minority_labels:
        minority_rows = train_df[train_df['v_call_light'].str.split('*').str[0] == label]
        upsampled_rows = minority_rows.sample(n=5000, replace=True)
        df_resampled = pd.concat([df_resampled, upsampled_rows])
    # Remove the minority rows from the original DataFrame
    df_majority = train_df[~train_df['v_call_light'].str.split('*').str[0].isin(minority_labels)]

    # Append the upsampled DataFrame to the original DataFrame
    df_resampled = pd.concat([df_majority, df_resampled])

    # Reset the index of the new DataFrame
    train_df = df_resampled.reset_index(drop=True)
    
    # create train and val set
    train_df['dataset'] = 'train set'
    val_df['dataset'] = 'val set'
    train_df.to_csv(f'{out_path}_train.tsv',sep='\t')
    val_df.to_csv(f'{out_path}_val.tsv',sep='\t')
    
    # create test set
    
    for i,t in enumerate(testset_ls):
        t['dataset'] = f'test set by {clusters[i]}'
    test_df = pd.concat(testset_ls)    
    test_df.to_csv(f'{out_path}_test.tsv',sep='\t')
    
    print('Training set: %d; Validation set: %d; Test set: %d' % (len(train_df),len(val_df),len(test_df)))




if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='clean dataset and split into train/test set')
    parser.add_argument('-i', '--input_file', default='results/OAS_memory_paired_v4.csv', type=str,
                        help='the input csv file')
    parser.add_argument('-o', '--output_prefix', default='data/dataset/OAS_memory', type=str,
                        help='the output xlsx file')
    parser.add_argument('-ci', '--cluster_identities',  default=['0.5','0.6','0.7','0.8','0.9'], nargs='+', help='cluster identity list to add the xlsx table')
    # Parse the arguments
    args = parser.parse_args()
    DF = pd.read_csv(args.input_file, low_memory=False)
    # output the df
    data_split(DF,args.cluster_identities,args.output_prefix)

