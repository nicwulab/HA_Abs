#!/usr/bin/python
import os
import sys
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict, Counter

def reading_cluster_info(filename):
  df = pd.read_csv(filename, sep="\t")
  cluster_dict = defaultdict(list)
  for index, info_dict in df.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    cluster_ID = info_dict['cluster']
    CDRH3      = info_dict['CDRH3_AA']
    cluster_dict[cluster_ID].append(CDRH3)
  return cluster_dict

def make_sequence_logo(sequence_list, figname):
  CDRH3_len = len(sequence_list[0])
  logo_width = (CDRH3_len-0)*0.6
  logo_width = 11 if logo_width > 11 else logo_width
  fig, ax = plt.subplots(1,1,figsize=[logo_width,2])
  seqlogo_matrix = logomaker.alignment_to_matrix(sequence_list)
  seqlogo = logomaker.Logo(seqlogo_matrix, font_name="Arial", color_scheme="weblogo_protein", width=1, ax=ax)
  seqlogo.style_spines(visible=False)
  seqlogo.ax.set_xticks([])
  seqlogo.ax.set_yticks([])
  seqlogo.fig.tight_layout()
  plt.savefig(figname, dpi=150)
  plt.close()
  print('Written %s' % figname, file = sys.stdout)

def write_CDRH3_seqlogo(cluster_dict, out_folder):
  for cluster_ID in sorted(cluster_dict.keys(), key=lambda x:int(x)):
    file_seqlogo = out_folder+'/seqlogo_'+cluster_ID+'.png'
    CDRH3s    = cluster_dict[cluster_ID]
    make_sequence_logo(CDRH3s, file_seqlogo)
    
def main():
  filename   = 'result/Ab_info_CDRH3_clustering.tsv'
  out_folder = 'CDRH3_seqlogo'
  os.system('mkdir -p '+out_folder)
  cluster_dict = reading_cluster_info(filename)
  write_CDRH3_seqlogo(cluster_dict, out_folder)

if __name__ == "__main__":
  main()
