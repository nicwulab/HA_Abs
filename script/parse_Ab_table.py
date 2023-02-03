#!/usr/bin/python
import os
import sys
import glob
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict, Counter

def extract_basic_info(filename, outfile_CDRH3, outfile_ref):
  aas = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  out_tsv = 'data/'+filename.rsplit('/')[1].replace('.xlsx','.tsv')
  print ('reading: %s' % filename)
  print ('writing: %s' % outfile_CDRH3)
  print ('writing: %s' % outfile_ref)
  print ('writing: %s' % out_tsv)
  outfile1 = open(outfile_CDRH3, 'w')
  datasheet = pd.read_excel(filename, sheet_name='Sheet1')
  datasheet.to_csv(out_tsv, sep="\t", index=False)
  Abs  = []
  refs = []
  donors = []
  epi    = []
  count_flu_Ab = 0
  for index, info_dict in datasheet.iterrows():
    for i in info_dict.keys():
      info_dict[i] = str(info_dict[i])
    Ab  = info_dict['Name'].replace('-','_')
    ref = info_dict['Reference']
    if Ab in Abs:
      print ("Name duplication: %s" % Ab)
    Abs.append(Ab)
    refs += [ref for ref in ref.rsplit('; ') if 'Unpublished' not in ref]
    donors.append(info_dict['Donor ID'])
    count_flu_Ab += 1
    epi.append(str(info_dict['Binds to'])) 
    if str(info_dict['Donor ID']) == 'nan': continue
    if str(info_dict['CDRH3_AA']) == 'nan': continue
    ID     = str(info_dict['Donor ID'])+'|'+str(info_dict['Name'])
    CDRH3  = str(info_dict['CDRH3_AA'])
    CDRL3  = str(info_dict['CDRL3_AA'])
    if 'PDB' in ref: print (ID)
    if '*' in CDRH3 or 'NA' == CDRH3: continue
    outfile1.write(ID+"\t"+CDRH3+"\n")
  outfile1.close()
  outfile3 = open(outfile_ref, 'w')
  outfile3.write("\n".join(sorted(list(set((refs)))))+"\n")
  outfile3.close()
  print ('# of Research Papers: %i' % len(set([ref for ref in refs if 'Patent' not in ref])))
  print ('# of Patents: %i' % len(set([ref for ref in refs if 'Patent' in ref])))
  print ('# of total donors: %i' % len(set(donors)))
  print ('# of flu Abs: %i' % count_flu_Ab)
  print ('  ',Counter(epi))

def main():
  filename = 'doc/HA_Abs_v14.10.xlsx'
  outfile_CDRH3  = 'result/CDRH3.tsv'
  outfile_ref    = 'result/refs.txt'
  extract_basic_info(filename, outfile_CDRH3, outfile_ref)

if __name__ == "__main__":
  main()
