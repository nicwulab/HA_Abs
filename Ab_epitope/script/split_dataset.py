import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
from sklearn.model_selection import GroupShuffleSplit,train_test_split
from sklearn.utils import shuffle
'''
firstly, we remove "unknown" group sequence that is similar to other annotated groups (sequence identity 0.9)
#consider to do: adding 99% sequence identity from "Unknown" group as "Known"
secondly, get dataset with ["HA:Head", "HA:Stem","HIV", "S:NTD", "S:RBD", "S:S2", "Others"]
split into train/test/val by cluster ID 0.8
'''

def df2fasta(df,output_file,name,seq):
    '''
    write a dataframe to fasta file
    :param df:
    :param output_file:
    :param name:
    :param seq:
    :return:
    '''
    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Loop over the rows in the DataFrame
        for index, row in df.iterrows():
            # Write the header line with the sequence name
            f.write('>{}\n'.format(row[name]))
            # Write the sequence to the file
            f.write('{}\n'.format(row[seq]))

# turn off the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

def data_clean(input):
    '''
    clean the epitope info table by cleaning and down-sampling Unknown category and up sampling other functional ones
    :param input: path to epitope info table
    :return: clean DF
    '''
    # data cleaning and data up/down-sampling
    fun_ls = ["HA:Head", "HA:Stem","HIV", "S:NTD", "S:RBD", "S:S2"]
    up_ls = ["HA:Head", "HA:Stem", "S:NTD", "S:S2"]
    df = pd.read_excel(input)
    clstr = list(set(df['cluster_ID_by_0.9']))
    dfs=[]
    print("start data cleaning \n")
    for c in tqdm(clstr):
        i_df = df[df['cluster_ID_by_0.9']==c]
        ## to remove "unknown" sequence similar to known
        # condition on same cluster, "Unknown" co-exist in other functions, label "Unknown" as T and drop them
        if ('Unknown' in set(i_df['Antigen_epitopes'])) & any(epi in fun_ls for epi in set(i_df['Antigen_epitopes'])):
            i_df.loc[:,'remove_unk'] = np.where(i_df['Antigen_epitopes']=='Unknown', 'T', 'F')
            # remove column `remove_unk` == 'T'
            i_df = i_df.drop(i_df[i_df['remove_unk'] == 'T'].index)
        # if all are "Unknown", randomly pick one
        elif ('Unknown' in set(i_df['Antigen_epitopes'])) & (len(set(i_df['Antigen_epitopes']))==1):
            i_df = i_df.sample(n=1)
        # if those "functions" are needed up sampling, randomly up-sampling
        elif any(epi in up_ls for epi in set(i_df['Antigen_epitopes'])):
            # set up sampling fold
            fold=20
            i_df = i_df.sample(n=fold*len(i_df),replace=True)

        # remove column 'Antigen_epitopes' == "HA:Unk"
        i_df = i_df.drop(i_df[i_df['Antigen_epitopes'] == 'HA:Unk'].index)
        dfs.append(i_df)
    DF= pd.concat(dfs)
    # replace "Unknown" with "Others"
    DF['Antigen_epitopes'] = DF.Antigen_epitopes.replace(to_replace='Unknown', value='Others')
    
    return DF
    
    

def data_split(df,cluster,out_path):
    '''
    split data into train/val/test by sequence identity 80%
    :param df:
    :param out_path:
    :return:
    '''
    print("\n start data splitting \n")
    df = df.sample(frac=1).reset_index(drop=True)
    splitter = GroupShuffleSplit(test_size=.10, n_splits=2, random_state=42)
    split = splitter.split(df, groups=df[cluster])
    train_inds, test_inds = next(split)
    train_df = df.iloc[train_inds]
    test_df = df.iloc[test_inds]

    train_df, val_df = train_test_split(train_df, test_size=0.1)
    val_df = val_df.drop_duplicates(keep='last')
    test_df = test_df.drop_duplicates(keep='last')
    train_df['dataset'] = 'train set'
    val_df['dataset'] = 'val set'
    test_df['dataset'] = 'test set'
    print('Training set: %d; Validation set: %d; Test set: %d' % (len(train_df),len(val_df),len(test_df)))
    train_df.to_csv(f'{out_path}_train.tsv',sep='\t')
    val_df.to_csv(f'{out_path}_val.tsv',sep='\t')
    test_df.to_csv(f'{out_path}_test.tsv',sep='\t')


if __name__ == "__main__":
    # Set up argparse
    parser = argparse.ArgumentParser(description='clean dataset and split into train/test set')
    parser.add_argument('-i', '--input_file', default='result/epitope_info_clstr_v2.xlsx', type=str,
                        help='the input Excel file')
    parser.add_argument('-o', '--output_prefix', default='result/epitope', type=str,
                        help='the output xlsx file')

    parser.add_argument('-c', '--clstr_col', default='cluster_ID_by_0.8', type=str, help="choose cluster identity column to split dataset")
    # Parse the arguments
    args = parser.parse_args()
    DF = data_clean(args.input_file)
    # output the df and fasta
    DF.to_excel(f'{args.output_prefix}_for_training_v1.xlsx',index=False)
    df2fasta(DF,f'{args.output_prefix}_clean.fasta','Name','VH_AA')
    data_split(DF,args.clstr_col,args.output_prefix)

