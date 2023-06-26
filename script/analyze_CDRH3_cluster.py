#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from collections import defaultdict, Counter

def reading_cluster_info(file_cluster):
  infile = open(file_cluster, 'r')
  cluster_dict = defaultdict(list)
  for line in infile.readlines():
    if 'cluster' in line: continue
    ID, CDRH3, cluster = line.rstrip().rsplit("\t")
    cluster_dict[cluster].append(ID)
  infile.close()
  return cluster_dict

def cleanup_cluster(cluster_dict):
  for n in sorted(cluster_dict.keys()):
    donors = list(set([ab.rsplit('|')[0] for ab in cluster_dict[n]]))
    if len(donors) == 1:
      del cluster_dict[n]
    elif len(donors) == 2 and 'VRC_310_cohort' in donors and ''.join(donors).count('VRC_310_cohort') == 2:
      del cluster_dict[n]
    elif len(donors) == 2 and 'Harrison_donors' in donors and ''.join(donors).count('Harrison_donor') == 2:
      del cluster_dict[n]
  return cluster_dict

def reading_Ab_table(filename):
  print ('reading: %s' % filename)
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  Ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    if str(info_dict['Donor ID']) == '': continue
    if str(info_dict['CDRH3_AA']) == '': continue
    ID = str(info_dict['Donor ID'])+'|'+str(info_dict['Name'])
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    Ab_dict[ID] = info_dict
  return Ab_dict, list(datasheet.columns)

def classify_epitope(epitopes):
  epi_list = []
  for epitope in epitopes:
    if epitope in ['NP', 'NA:Unk', 'Other'] or epitope == '':
      continue
    elif epitope in ['HA:Unk', 'H1 HA', 'Multiple', 'Unknown', 'HA or NP or M1']:
      epi_list.append('HA:Stem')
      epi_list.append('HA:Head')
    elif epitope in ['HA:Stem', 'HA:Stem, CR9114-competing', 'HA:Stem, non-CR9114-competing']:
      epi_list.append('HA:Stem')
    elif epitope == 'HA:Head':
      epi_list.append('HA:Head')
    else:
      print (epitope)
  if len(epi_list) == 0:
    return ('unknown')
  max_epi = max(epi_list,key=epi_list.count)
  if epi_list.count(max_epi)/float(len(epi_list)) > 0.5:
    return (max_epi)
  else:
    return ('unknown')

def merge_Vgene(Vgene):
  if 'IGHV4-38' in Vgene: 
    return (Vgene)
  if 'IGHV5-10' in Vgene: 
    return (Vgene)
  if Vgene == 'IGHV3-53' or Vgene == 'IGHV3-66':
    return 'IGHV3-53/3-66'
  if Vgene == 'IGHV4-30-2' or Vgene == 'IGHV4-30-4':
    return 'IGHV4-30-2/4-30-4'
  if Vgene == 'IGHV3-30-3':
    return 'IGHV3-30'
  if Vgene.count('-') == 2:
    adj_Vgene = '-'.join(Vgene.rsplit('-')[0:2])
    print ('converting: %s to %s' % (Vgene, adj_Vgene))
    return adj_Vgene
  else:
    return Vgene

def find_top_Vgene(Vgenes):
  Vgenes = list(map(merge_Vgene, [Vgene.rsplit(',')[0].rsplit('*')[0] for Vgene in Vgenes]))
  top_Vgene = max(Vgenes,key=Vgenes.count)
  freq      = Vgenes.count(top_Vgene)/float(len(Vgenes))
  if len(set(Vgenes))==1 and Vgenes[0]=='nan':
    return 'nan', 'nan'
  elif freq <= 0.5:
    return 'Diverse', freq
  else:
    return top_Vgene, freq

def extract_Vgene_info(HC_Vgenes, LC_Vgenes):
  top_HV, top_HV_freq = find_top_Vgene(HC_Vgenes)
  top_LV, top_LV_freq = find_top_Vgene(LC_Vgenes)
  return top_HV, top_HV_freq, top_LV, top_LV_freq

def seqs2consensus(seqlist):
  consensus = ''
  for n in range(len(seqlist[0])):
    resi = []
    for seq in seqlist:
      resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

def write_cluster_info(outfile1, outfile2, header, Ab_dict, cluster_dict):
  print ('writing: %s' % outfile1)
  print ('writing: %s' % outfile2)
  outfile1 = open(outfile1, 'w')
  outfile2 = open(outfile2, 'w')
  outfile1.write('cluster'+"\t"+"\t".join(header)+"\n")
  outfile2.write("\t".join(['cluster_ID', 'cluster_size', 'num_donors', 'epitope', 
                            'top_HV', 'top_LV', 'top_HV_freq', 'top_LV_freq', 'CDRH3_consen', 'donors', 'ab_names'])+"\n")
  cluster_ID = 0
  for n in sorted(cluster_dict.keys(), key=lambda x:len(cluster_dict[x]), reverse=True):
    if str(Ab_dict[cluster_dict[n][0]]['CDRH3_AA']) in ['NA', 'nan']:
      print ("One cluster with 'NA' CDRH3 is removed")
      continue
    cluster_ID += 1
    HC_Vgenes = []
    LC_Vgenes = []
    CDRH3s    = []
    donor_names  = []
    ab_names  = []
    for Ab in cluster_dict[n]:
      outfile1.write(str(cluster_ID)+"\t"+"\t".join([str(Ab_dict[Ab][str(item)]) for item in header])+"\n")
      HC_Vgenes.append(str(Ab_dict[Ab]['Heavy_V_gene']))
      LC_Vgenes.append(str(Ab_dict[Ab]['Light_V_gene']))
      CDRH3s.append(str(Ab_dict[Ab]['CDRH3_AA']))
      donor_names.append(Ab.rsplit('|')[0])
      ab_names.append(Ab.rsplit('|')[1])
    CDRH3_consen = seqs2consensus(CDRH3s)
    top_HV, top_HV_freq, top_LV, top_LV_freq = extract_Vgene_info(HC_Vgenes, LC_Vgenes)
    cluster_size = len(cluster_dict[n])
    donors       = len(set([Ab.rsplit('|')[0] for Ab in cluster_dict[n]]))
    epitopes     = [Ab_dict[Ab]['Binds to'] for Ab in cluster_dict[n]]
    epitope      = classify_epitope(epitopes)
    outfile2.write("\t".join(map(str, [cluster_ID, cluster_size, donors, epitope,
                                       top_HV, top_LV, top_HV_freq, top_LV_freq, CDRH3_consen,
                                       ','.join(sorted(list(set(donor_names)))),
                                       ','.join(sorted(list(set(ab_names))))]))+"\n")
  outfile1.close()
  outfile2.close()

def main():
  outfile1  = 'result/Ab_info_CDRH3_clustering.tsv'
  outfile2  = 'result/CDRH3_cluster_summary.tsv'
  filename = 'doc/HA_Abs_v17.xlsx'
  file_cluster  = 'result/CDRH3_cluster.tsv'
  Ab_dict, header  = reading_Ab_table(filename)
  cluster_dict  = reading_cluster_info(file_cluster)
  print ("total number of clusters: %i" % len(cluster_dict.keys()))
  cluster_dict  = cleanup_cluster(cluster_dict)
  print ("total number of clusters after clean up: %i" % len(cluster_dict.keys()))
  write_cluster_info(outfile1, outfile2, header, Ab_dict, cluster_dict)

if __name__ == "__main__":
  main()
