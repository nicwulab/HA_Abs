#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
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

def extract_ab(Ab_dict, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['name', 'donor', 'HV', 'LV', 'CDRH3', 'YGD', 'PDB', "ref"])+"\n")
  for Ab in Ab_dict.keys():
    donor  = Ab_dict[Ab]['Donor ID']
    epi    = Ab_dict[Ab]['Binds to']
    PDB    = Ab_dict[Ab]['PDB']
    VH_nuc = Ab_dict[Ab]['VH_nuc']
    HV     = Ab_dict[Ab]['Heavy_V_gene'].rsplit('*')[0]
    LV     = Ab_dict[Ab]['Light_V_gene'].rsplit('*')[0]
    HDs    = list(set([HD.rsplit('*')[0] for HD in Ab_dict[Ab]['Heavy_D_gene'].rsplit(',')]))
    CDRH3  = Ab_dict[Ab]['CDRH3_AA']
    ref    = Ab_dict[Ab]['Reference']
    if epi == 'HA:Head' and len(HDs) == 1 and HDs[0] == 'IGHD4-17':
      if 'YGD' in CDRH3:
        YGD = 'yes'
      elif any(x in CDRH3 for x in ['YAD', 'YGE', 'FGD']):
        YGD = 'similar'
      else: 
        YGD = 'no'
      outfile.write("\t".join([Ab, donor, HV, LV, CDRH3, YGD, PDB, ref])+"\n")
  outfile.close()

def main():
  filename = 'doc/HA_Abs_v18.xlsx'
  outfile  = 'result/Abs_IGHD4-17.tsv'
  Ab_dict  = read_ab_info(filename)
  extract_ab(Ab_dict, outfile)

if __name__ == "__main__":
  main()
