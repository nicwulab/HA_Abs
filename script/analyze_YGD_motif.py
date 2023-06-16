#!/usr/bin/python
import os
import sys
import glob
import math
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

def motif_freq(CDRH3s, motif):
  total_count = float(len(CDRH3s))
  freq = float(len([CDRH3 for CDRH3 in CDRH3s if motif in CDRH3]))/total_count
  SE = math.sqrt(freq*(1-freq)/(total_count))
  return {'freq':freq, 'SE':SE}

def write_motif_freq_file(motif_freq_dict, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['Ab_type', 'freq', 'SE'])+"\n")
  for ab_type in motif_freq_dict.keys():
    freq = motif_freq_dict[ab_type]['freq']
    SE   = motif_freq_dict[ab_type]['SE']
    outfile.write("\t".join(map(str,[ab_type, freq, SE]))+"\n")
  outfile.close()

def main():
  outfile  = 'result/YGD_motif_freq.tsv'
  filename = 'doc/HA_Abs_v17.xlsx'
  file_GB  = 'doc/all_paired_antibodies_from_GB_v6.xlsx'
  Ab_dict  = read_ab_info(filename)
  GB_dict  = read_ab_info(file_GB)
  motif    = 'YGD'
  stem_CDRH3s        = [Ab_dict[Ab]['CDRH3_AA'] for Ab in Ab_dict.keys() if Ab_dict[Ab]['Binds to'] == 'HA:Stem']
  head_CDRH3s        = [Ab_dict[Ab]['CDRH3_AA'] for Ab in Ab_dict.keys() if Ab_dict[Ab]['Binds to'] == 'HA:Head']
  head_non270_CDRH3s = [Ab_dict[Ab]['CDRH3_AA'] for Ab in Ab_dict.keys() \
                        if Ab_dict[Ab]['Binds to'] == 'HA:Head' and Ab_dict[Ab]['Heavy_V_gene'].rsplit('*')[0] != 'IGHV2-70']
  GB_CDRH3s          = [GB_dict[Ab]['CDRH3_AA'] for Ab in GB_dict.keys()]
  motif_freq_dict = {'Head':motif_freq(head_CDRH3s, motif),
                     'Head (non-IGHV2-70)':motif_freq(head_non270_CDRH3s, motif),
                     'Stem':motif_freq(stem_CDRH3s, motif),
                     'GenBank':motif_freq(GB_CDRH3s, motif),
                     }
  write_motif_freq_file(motif_freq_dict, outfile)
  print ('Total # of CDRH3 for HA stem: %i' % len(stem_CDRH3s))
  print ('Total # of CDRH3 for HA head: %i' % len(head_CDRH3s))
  print ('Total # of CDRH3 for HA head (non-IGHV2-70): %i' % len(head_non270_CDRH3s))
  print ('Total # of CDRH3 for GenBank: %i' % len(GB_CDRH3s))
  
if __name__ == "__main__":
  main()
