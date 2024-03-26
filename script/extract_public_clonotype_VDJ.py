#!/usr/bin/python
import os
import sys
import glob
import math
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict, Counter

def read_ab_info(filename):
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

def identifying_max_gene(gene_list):
  return max(gene_list,key=gene_list.count)

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

def IG_list_to_freq(gene_dict):
  for epi in gene_dict.keys():
    gene_dict[epi] = Counter(gene_dict[epi])
    if 'nan' in gene_dict[epi].keys():
      del gene_dict[epi]['nan']
    total_count = float(sum(gene_dict[epi].values()))
    for IG in gene_dict[epi].keys():
      freq = float(gene_dict[epi][IG])/total_count
      SE = math.sqrt(freq*(1-freq)/(total_count))
      gene_dict[epi][IG] = {'freq':round(freq,4), 'SE':round(SE,4)}
  return gene_dict

def IG_gene_analysis(Ab_dict, clono_dict):
  IGHV_gene_dict = {'HA:Stem':[], 'HA:Head':[]}
  IGHD_gene_dict = {'HA:Stem':[], 'HA:Head':[]}
  IGHJ_gene_dict = {'HA:Stem':[], 'HA:Head':[]}
  IGLV_gene_dict = {'HA:Stem':[], 'HA:Head':[]}
  IGLJ_gene_dict = {'HA:Stem':[], 'HA:Head':[]}
  IGV_pair_dict  = {'HA:Stem':[], 'HA:Head':[]}
  IGHV_clonotype     = defaultdict(list)
  IGHD_clonotype     = defaultdict(list)
  IGHJ_clonotype     = defaultdict(list)
  IGLV_clonotype     = defaultdict(list)
  IGLJ_clonotype     = defaultdict(list)
  IGV_pair_clonotype = defaultdict(list)
  epi_clonotype      = defaultdict(list)
  clono_IDs      = []
  for Ab in Ab_dict.keys():
    HV = Ab_dict[Ab]['Heavy_V_gene'].rsplit('*')[0]
    LV = Ab_dict[Ab]['Light_V_gene'].rsplit('*')[0]
    HDs = [HD.rsplit('*')[0] for HD in Ab_dict[Ab]['Heavy_D_gene'].rsplit(',')]
    HJs = [HJ.rsplit('*')[0] for HJ in Ab_dict[Ab]['Heavy_J_gene'].rsplit(',')]
    LJs = [LJ.rsplit('*')[0] for LJ in Ab_dict[Ab]['Light_J_gene'].rsplit(',')]
    donor = Ab_dict[Ab]['Donor ID']  
    epi   = Ab_dict[Ab]['Binds to']
    epi   = 'HA:Stem' if epi in ['HA:Stem, CR9114-competing', 'HA:Stem, non-CR9114-competing'] else epi
    clono = clono_dict[Ab] if Ab in clono_dict.keys() else 'none'
    clono_ID = clono+'|'+donor
    clono_IDs.append(clono_ID)
    if clono == 'none':
      if epi in ['HA:Stem', 'HA:Head']:
        IGHV_gene_dict[epi].append(HV)
        IGLV_gene_dict[epi].append(LV)
        if len(set(HDs)) == 1:
          IGHD_gene_dict[epi].append(HDs[0])
        if len(set(HJs)) == 1:
          IGHJ_gene_dict[epi].append(HJs[0])
        if len(set(LJs)) == 1:
          IGLJ_gene_dict[epi].append(LJs[0])
        if HV != 'nan' and LV != 'nan':
          IGV_pair_dict[epi].append(HV+'__'+LV)
    else:
      IGHV_clonotype[clono_ID].append(HV)
      IGLV_clonotype[clono_ID].append(LV)
      IGV_pair_clonotype[clono_ID].append(HV+'__'+LV)
      epi_clonotype[clono_ID].append(epi)
      if len(set(HDs)) == 1:
        IGHD_clonotype[clono_ID].append(HDs[0])
      if len(set(HJs)) == 1:
        IGHJ_clonotype[clono_ID].append(HJs[0])
      if len(set(LJs)) == 1:
        IGLJ_clonotype[clono_ID].append(LJs[0])
  for clono_ID in sorted(IGHV_clonotype.keys(), key=lambda x:float(x.rsplit('|')[0])):
    epi = classify_epitope(epi_clonotype[clono_ID.rsplit('|')[0]])
    if epi in ['HA:Stem', 'HA:Head']:
      IGHV = identifying_max_gene(IGHV_clonotype[clono_ID])
      IGLV = identifying_max_gene(IGLV_clonotype[clono_ID])
      V_pair = identifying_max_gene(IGV_pair_clonotype[clono_ID])
      IGHV_gene_dict[epi].append(IGHV)
      IGLV_gene_dict[epi].append(IGLV)
      IGV_pair_dict.append(V_pair)
      if len(IGHD_clonotype[clono_ID]) != 0:
        IGHD = identifying_max_gene(IGHD_clonotype[clono_ID])
        IGHD_gene_dict[epi].append(IGHD)
      if len(IGHJ_clonotype[clono_ID]) != 0:
        IGHJ = identifying_max_gene(IGHJ_clonotype[clono_ID])
        IGHJ_gene_dict[epi].append(IGHJ)
      if len(IGLJ_clonotype[clono_ID]) != 0:
        IGLJ = identifying_max_gene(IGLJ_clonotype[clono_ID])
        IGLJ_gene_dict[epi].append(IGLJ)
  return sorted(list(set(clono_IDs))), IGHV_gene_dict, IGHD_gene_dict, IGLV_gene_dict, IGHJ_gene_dict, IGLJ_gene_dict, IGV_pair_dict
 
def VDJ_baseline(GB_dict, IGHV_gene_dict, IGHD_gene_dict, IGHJ_gene_dict, IGLV_gene_dict, IGLJ_gene_dict, IGV_pair_dict):
  IGHV_gene_dict['GenBank'] = []
  IGHD_gene_dict['GenBank'] = []
  IGHJ_gene_dict['GenBank'] = []
  IGLV_gene_dict['GenBank'] = []
  IGLJ_gene_dict['GenBank'] = []
  IGV_pair_dict['GenBank']  = []
  for Ab in GB_dict.keys():
    HV = GB_dict[Ab]['Heavy_V_gene'].rsplit('*')[0]
    LV = GB_dict[Ab]['Light_V_gene'].rsplit('*')[0]
    IGHV_gene_dict['GenBank'].append(HV)
    IGLV_gene_dict['GenBank'].append(LV)
    IGV_pair_dict['GenBank'].append(HV+'__'+LV)
    HDs = [HD.rsplit('*')[0] for HD in GB_dict[Ab]['Heavy_D_gene'].rsplit(',')]
    HJs = [HJ.rsplit('*')[0] for HJ in GB_dict[Ab]['Heavy_J_gene'].rsplit(',')]
    LJs = [LJ.rsplit('*')[0] for LJ in GB_dict[Ab]['Light_J_gene'].rsplit(',')]
    if len(set(HDs)) == 1:
      IGHD_gene_dict['GenBank'].append(HDs[0])
    if len(set(HJs)) == 1:
      IGHJ_gene_dict['GenBank'].append(HJs[0])
    if len(set(LJs)) == 1:
      IGLJ_gene_dict['GenBank'].append(LJs[0])
  return IGHV_gene_dict, IGHD_gene_dict, IGLV_gene_dict , IGHJ_gene_dict, IGLJ_gene_dict, IGV_pair_dict

def write_IG_freq(gene_dict, GB_dict, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['epitope', 'gene', 'freq', 'SE'])+"\n")
  gene_dict = IG_list_to_freq(gene_dict)
  IG_genes = list(set([IG_gene for epi in ['HA:Head','HA:Stem'] for IG_gene in gene_dict[epi].keys()]))
  for epi in sorted(gene_dict.keys()):
    for IG_gene in sorted(IG_genes):
      epi_rename = epi.rsplit('HA:')[1] if 'HA:' in epi else epi
      freq = gene_dict[epi][IG_gene]['freq'] if gene_dict[epi][IG_gene] != 0 else 0
      SE   = gene_dict[epi][IG_gene]['SE'] if gene_dict[epi][IG_gene] != 0 else 0
      outfile.write("\t".join(map(str, [epi_rename, IG_gene, freq, SE]))+"\n")
  outfile.close()

def identifying_public_clono(filename, outfile, clono_IDs):
  print ('writing: %s' % outfile)
  clono_count = Counter([clono_ID.rsplit('|')[0] for clono_ID in clono_IDs])
  clono_pub = [float(clono_ID) for clono_ID in clono_count.keys() if clono_count[clono_ID] > 1 and clono_ID != 'none']
  print (sorted(clono_pub))
  df = pd.read_excel(filename, sheet_name="Sheet1")
  df = df[df['clonotype'].isin(clono_pub)]
  print ('# of public clonotype %i' % len(clono_pub))
  df.to_excel(outfile, index=False)

def main():
  out_IGHV_file  = 'result/IGHV_freq.tsv'
  out_IGHD_file  = 'result/IGHD_freq.tsv'
  out_IGLV_file  = 'result/IGLV_freq.tsv'
  out_IGHJ_file  = 'result/IGHJ_freq.tsv'
  out_IGLJ_file  = 'result/IGLJ_freq.tsv'
  out_Vpair_file = 'result/IGV_pair_freq.tsv'
  out_pub_clono  = 'result/HA_Abs_clonotype_public.xlsx'
  file_ab        = 'doc/HA_Abs_v18.xlsx'
  file_clonotype = 'result/HA_Abs_clonotype.xlsx'
  file_GB        = 'doc/all_paired_antibodies_from_GB_v6.xlsx'
  Ab_dict        = read_ab_info(file_ab)
  GB_dict        = read_ab_info(file_GB)
  clono_dict     = read_clonoinfo(file_clonotype)
  clono_IDs, IGHV_gene_dict, IGHD_gene_dict, IGLV_gene_dict, IGHJ_gene_dict, IGLJ_gene_dict, IGV_pair_dict = IG_gene_analysis(Ab_dict, \
                                                                                                                              clono_dict)
  IGHV_gene_dict, IGHD_gene_dict, IGLV_gene_dict, IGHJ_gene_dict, IGLJ_gene_dict, IGV_pair_dict = VDJ_baseline(GB_dict, IGHV_gene_dict, \
                                                                                               IGHD_gene_dict, IGHJ_gene_dict, \
                                                                                               IGLV_gene_dict, IGLJ_gene_dict, \
                                                                                               IGV_pair_dict)
  write_IG_freq(IGHV_gene_dict, GB_dict, out_IGHV_file)
  write_IG_freq(IGHD_gene_dict, GB_dict, out_IGHD_file)
  write_IG_freq(IGLV_gene_dict, GB_dict, out_IGLV_file)
  write_IG_freq(IGHJ_gene_dict, GB_dict, out_IGHJ_file)
  write_IG_freq(IGLJ_gene_dict, GB_dict, out_IGLJ_file)
  write_IG_freq(IGV_pair_dict, GB_dict, out_Vpair_file)
  identifying_public_clono(file_clonotype, out_pub_clono, clono_IDs)

if __name__ == "__main__":
  main()
