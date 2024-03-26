#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats
from Bio.Seq import Seq
from collections import defaultdict, Counter

def read_ab_info(filename):
  aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  print ('reading: %s' % filename)
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  Ab_dict = {}
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    Ab_dict[info_dict['Name']] = info_dict
  return Ab_dict

def read_clonoinfo(filename):
  print ('reading: %s' % filename)
  df = pd.read_excel(filename, sheet_name="Sheet1")
  df = df[df.clonotype.notnull()]
  clono_dict = {}
  for index, info_dict in df.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    clono_dict[info_dict['Name']] = info_dict['clonotype']
  return clono_dict

def classify_epitope(epitopes):
  epi_list = []
  for epitope in epitopes:
    if epitope in ['NP', 'NA:Unk', 'Other'] or epitope == '':
      continue
    elif epitope in ['HA:Unk']:
      epi_list.append('HA:Stem')
      epi_list.append('HA:Head')
    elif epitope in ['HA:Stem']:
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

def seqs2consensus(seqlist, most_common_len):
  consensus = ''
  for n in range(most_common_len):
    resi = []
    for seq in seqlist:
      if len(seq) == most_common_len:
        resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

def generate_CDR_concensus(CDR_dict):
  for clono_ID in CDR_dict.keys():
    CDR_lens = [len(CDR) for CDR in CDR_dict[clono_ID]]
    most_common_len, num_most_common_len = Counter(CDR_lens).most_common(1)[0]
    CDR_dict[clono_ID] = seqs2consensus(CDR_dict[clono_ID], most_common_len)
  return CDR_dict

def compress_ab_list_CDRH3(Ab_dict, clono_dict):
  CDRH1_dict = defaultdict(list)
  CDRH2_dict = defaultdict(list)
  CDRH3_dict = defaultdict(list)
  CDRL1_dict = defaultdict(list)
  CDRL2_dict = defaultdict(list)
  CDRL3_dict = defaultdict(list)
  VH_dict    = defaultdict(list)
  VL_dict    = defaultdict(list)
  epi_dict = defaultdict(list)
  for Ab in Ab_dict.keys():
    CDRH1 = Ab_dict[Ab]['CDRH1_AA']
    CDRH2 = Ab_dict[Ab]['CDRH2_AA']
    CDRH3 = Ab_dict[Ab]['CDRH3_AA']
    CDRL1 = Ab_dict[Ab]['CDRL1_AA']
    CDRL2 = Ab_dict[Ab]['CDRL2_AA']
    CDRL3 = Ab_dict[Ab]['CDRL3_AA']
    VH    = Ab_dict[Ab]['VH_AA']
    VL    = Ab_dict[Ab]['VL_AA']
    donor = Ab_dict[Ab]['Donor ID']
    clono = clono_dict[Ab] if Ab in clono_dict.keys() else 'none'
    clono_ID = Ab if clono=='none' else clono+'|'+donor 
    epi      = Ab_dict[Ab]['Binds to']
    epi_dict[clono_ID].append(epi)
    if CDRH1 != 'nan': CDRH1_dict[clono_ID].append(CDRH1)
    if CDRH2 != 'nan': CDRH2_dict[clono_ID].append(CDRH2)
    if CDRH3 != 'nan': CDRH3_dict[clono_ID].append(CDRH3)
    if CDRL1 != 'nan': CDRL1_dict[clono_ID].append(CDRL1)
    if CDRL2 != 'nan': CDRL2_dict[clono_ID].append(CDRL2)
    if CDRL3 != 'nan': CDRL3_dict[clono_ID].append(CDRL3)
    if VH != 'nan': VH_dict[clono_ID].append(VH)
    if VL != 'nan': VL_dict[clono_ID].append(VL)
  for clono_ID in epi_dict.keys():
    epi_dict[clono_ID]  = classify_epitope(epi_dict[clono_ID])
  CDRH1_dict = generate_CDR_concensus(CDRH1_dict)
  CDRH2_dict = generate_CDR_concensus(CDRH2_dict)
  CDRH3_dict = generate_CDR_concensus(CDRH3_dict)
  CDRL1_dict = generate_CDR_concensus(CDRL1_dict)
  CDRL2_dict = generate_CDR_concensus(CDRL2_dict)
  CDRL3_dict = generate_CDR_concensus(CDRL3_dict)
  VH_dict = generate_CDR_concensus(VH_dict)
  VL_dict = generate_CDR_concensus(VL_dict)
  return CDRH1_dict, CDRH2_dict, CDRH3_dict, CDRL1_dict, CDRL2_dict, CDRL3_dict, VH_dict, VL_dict, epi_dict

def add_GB_ab(CDRH3_dict, epi_dict, GB_dict):
  for Ab in GB_dict.keys():
    CDRH3 = GB_dict[Ab]['CDRH3_AA']
    if CDRH3 != 'nan': CDRH3_dict[Ab] = CDRH3
    epi_dict[Ab] = 'GenBank'
  return CDRH3_dict, epi_dict

def analyze_CDRH3_property(CDRH3_dict, epi_dict, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['clone', 'epi', 'CDRH3', 'CDRH3_tip', 'len_CDRH3',
                           'CDRH3_H_score', 'CDRH3_C_score', 'CDRH3_A_score',
                           'tip_H_score', 'tip_C_score', 'tip_A_score'])+"\n")
  H_scale = {'A':0.17,'R':0.81,'N':0.42,'D':1.23,'C':-0.24,
             'Q':0.58,'E':2.02,'G':0.01,'H':0.96,'I':-0.31,
             'L':-0.56,'K':0.99,'M':-0.23,'F':-1.13,'P':0.45,
             'S':0.13,'T':0.14,'W':-1.85,'Y':-0.94,'V':0.07}
  C_scale = {'A':0,'R':1,'N':0,'D':-1,'C':0,
             'Q':0,'E':-1,'G':0,'H':0.25,'I':0,
             'L':0,'K':1,'M':0,'F':0,'P':0,
             'S':0,'T':0,'W':0,'Y':0,'V':0}
  A_scale = {'A':0,'R':0,'N':0,'D':0,'C':0,
             'Q':0,'E':0,'G':0,'H':0,'I':0,
             'L':0,'K':0,'M':0,'F':1,'P':0,
             'S':0,'T':0,'W':1,'Y':1,'V':0}
  for clono_ID in CDRH3_dict.keys():
    CDRH3 = CDRH3_dict[clono_ID]
    epi   = epi_dict[clono_ID].replace('HA:','')
    if set(list(CDRH3)).issubset(H_scale.keys()):
      if len(CDRH3)%2==0:
        CDRH3_tip = CDRH3[int(len(CDRH3)/2)-2:int(len(CDRH3)/2)+2]
      else:
        CDRH3_tip = CDRH3[int(len(CDRH3)/2)-1:int(len(CDRH3)/2)+2]
      tip_H_score = -sum([H_scale[aa] for aa in CDRH3_tip])/float(len(CDRH3_tip))*10
      tip_C_score = sum([abs(C_scale[aa]) for aa in CDRH3_tip])/float(len(CDRH3_tip))
      tip_A_score = sum([A_scale[aa] for aa in CDRH3_tip])/float(len(CDRH3_tip))
      CDRH3_H_score = -sum([H_scale[aa] for aa in CDRH3])/float(len(CDRH3))*10
      CDRH3_C_score = sum([abs(C_scale[aa]) for aa in CDRH3])/float(len(CDRH3))
      CDRH3_A_score = sum([A_scale[aa] for aa in CDRH3])/float(len(CDRH3))
    outfile.write("\t".join(map(str, [clono_ID, epi, CDRH3, CDRH3_tip, len(CDRH3),
                  CDRH3_H_score, CDRH3_C_score, CDRH3_A_score,
                  tip_H_score, tip_C_score, tip_A_score]))+"\n")
  outfile.close()

def write_CDRs(CDRH1_dict, CDRH2_dict, CDRH3_dict, CDRL1_dict, CDRL2_dict, CDRL3_dict, VH_dict, VL_dict, epi_dict, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  clones = set(CDRH1_dict.keys()).intersection(set(CDRH2_dict.keys())).intersection(set(CDRH3_dict.keys())) \
           .intersection(set(CDRL1_dict.keys())).intersection(set(CDRL2_dict.keys())).intersection(set(CDRL3_dict.keys())) \
           .intersection(set(VH_dict)).intersection(VL_dict)
  outfile.write("\t".join(['Name', 'Binds to', 'VH_AA', 'VL_AA',
                           'CDRH1_AA', 'CDRH2_AA', 'CDRH3_AA', 'CDRL1_AA', 'CDRL2_AA', 'CDRL3_AA'])+"\n")
  for clone in sorted(clones, key=lambda x:epi_dict[x]):
    epi   = epi_dict[clone]
    VH_AA = VH_dict[clone]
    VL_AA = VL_dict[clone]
    CDRH1 = CDRH1_dict[clone]
    CDRH2 = CDRH2_dict[clone]
    CDRH3 = CDRH3_dict[clone]
    CDRL1 = CDRL1_dict[clone]
    CDRL2 = CDRL2_dict[clone]
    CDRL3 = CDRL3_dict[clone]
    if epi in ['HA:Stem', 'HA:Head']:
      outfile.write("\t".join([clone, epi, VH_AA, VL_AA, CDRH1, CDRH2, CDRH3, CDRL1, CDRL2, CDRL3])+"\n")
  outfile.close()

def main():
  filename       = 'doc/HA_Abs_v18.xlsx'
  outfile1       = 'result/CDRH3_property.tsv'
  outfile2       = 'result/Ab_for_model.tsv'
  file_clonotype = 'result/HA_Abs_clonotype.xlsx'
  file_GB        = 'doc/all_paired_antibodies_from_GB_v6.xlsx'
  Ab_dict        = read_ab_info(filename)
  GB_dict        = read_ab_info(file_GB)
  clono_dict     = read_clonoinfo(file_clonotype)
  CDRH1_dict,CDRH2_dict,CDRH3_dict,CDRL1_dict,CDRL2_dict,CDRL3_dict,VH_dict,VL_dict,epi_dict = compress_ab_list_CDRH3(Ab_dict, clono_dict)
  CDRH3_dict, epi_dict = add_GB_ab(CDRH3_dict, epi_dict, GB_dict)
  analyze_CDRH3_property(CDRH3_dict, epi_dict, outfile1)
  write_CDRs(CDRH1_dict, CDRH2_dict, CDRH3_dict, CDRL1_dict, CDRL2_dict, CDRL3_dict, VH_dict, VL_dict, epi_dict, outfile2)

if __name__ == "__main__":
  main()
