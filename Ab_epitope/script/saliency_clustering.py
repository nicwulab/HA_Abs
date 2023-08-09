import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
# load all dataset
DF = pd.read_excel('../result/epitope_info_clstr_v2.xlsx').drop_duplicates()
stem_df = DF[DF['Antigen_epitopes'] == 'HA:Stem']


gradcam_dict = defaultdict(list)
seq_dict = {}
gene_dict = {}
for idx,row in stem_df.iterrows():
    name = row['Name'].replace('/', '_')
    epi = row['Antigen_epitopes']
    file = f'{dir}{name}_gradcam.tsv'
    try:
        seq_dict[name]=row['VH_AA']
        gene_dict[name]=row['Heavy_V_gene'].split('*')[0]
        # Create an empty list to store each DataFrame
        dfs = []
        
        # Loop through each file and read it into a DataFrame
        for i in range(15):
            file_path = f'../dataset/explain_mBLM_{i+1}/{name}_gradcam.tsv'
            df = pd.read_csv(file_path, sep='\t').rename(columns={'0': name}).loc[:,name]
            dfs.append(df)
        # Concatenate all DataFrames into a single DataFrame
        combined_df = pd.concat(dfs,axis=1)
        # Calculate the mean value across the rows and store it in a new DataFrame
        df_saliency = combined_df.mean(axis=1) 
        df_saliency = df_saliency.rename(name)
        df_saliency.to_csv(f'../result/average_saliency/{name}_gradcam.tsv', sep='\t')
        gradcam_dict[epi].append(df_saliency)
    except FileNotFoundError:
        print(f"{file} not found. Skipping...")
        continue
        
from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as sch
# clean df
stem_concat = pd.concat(gradcam_dict['HA:Stem'], axis=1)
stem_concat = stem_concat.loc[:,~stem_concat.columns.duplicated()]
# Calculate the pairwise distances between rows
distance_matrix = hierarchy.distance.pdist(stem_concat.T.values)

# Perform hierarchical clustering
linkage_matrix = hierarchy.linkage(distance_matrix, method='ward', metric='euclidean')

# Obtain column cluster labels
row_clusters = sch.fcluster(linkage_matrix, t=6, criterion='maxclust')

# Create a list of DataFrames for each cluster
df = stem_concat.T
df['cluster']=row_clusters

sampled_df = df.groupby('cluster').apply(lambda x: x.sample(n=16, random_state=42))
# Drop column 'B'
sampled_df = sampled_df.drop('cluster', axis=1)


# Perform hierarchical clustering
col_linkage = sch.linkage(sampled_df.values, method='ward', metric='euclidean')


# make a cluster heatmap plot
p = sns.clustermap(sampled_df.T, cmap="coolwarm", col_linkage=col_linkage,
                   standard_scale=1, row_cluster=False,figsize=(3, 6),
                   tree_kws=dict(linewidths=1.5, colors='grey'))

p.ax_heatmap.set_ylim(120, 0)
# p.show
p.savefig(f'../graph/mBLM_saliency_clusters.png', dpi=300)