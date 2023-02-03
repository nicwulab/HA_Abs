#/usr/bin/python
import sys
import sklearn
import distance
from Bio import SeqIO
from collections import defaultdict

def file_to_dict(filename):
  CDRH3_dict = {}
  print ("reading: %s" % filename)
  infile = open(filename, 'r') 
  for line in infile.readlines():
    ID, CDRH3 = line.rstrip().rsplit("\t")
    CDRH3_dict[ID] = CDRH3
  return (CDRH3_dict)

def clustering_CDRH3(CDRH3_dict, cutoff):
  ID_list = list(sorted(CDRH3_dict.keys(),key=lambda x:CDRH3_dict[x]))
  clusters = {} #cluster[length][clusterID] = list of Ab name
  count_ID = 0
  for ID_a in ID_list:
    count_ID += 1
    if count_ID%100 == 0: print ('processed: %i Abs' % count_ID)
    CDRH3_a = CDRH3_dict[ID_a]
    if CDRH3_a == 'NA': continue
    CDRH3_a_len = len(CDRH3_a)
    assigned = 0
    if CDRH3_a_len in clusters.keys():
      for n in sorted(clusters[CDRH3_a_len].keys()):
        ID_b = clusters[CDRH3_a_len][n][0]
        CDRH3_b = CDRH3_dict[ID_b]
        dist = distance.hamming(CDRH3_a, CDRH3_b)
        if float(dist)/float(CDRH3_a_len) < cutoff:
          assigned = 1
          clusters[CDRH3_a_len][n].append(ID_a)
          break
    else:
      clusters[CDRH3_a_len] = defaultdict(dict)
    if assigned == 0:
      new_cluster_ID = len(clusters[CDRH3_a_len].keys())+1
      clusters[CDRH3_a_len][new_cluster_ID] = [ID_a]
  return (clusters)

def merging_cluster(clusters, CDRH3_dict, cutoff):
  for CDRH3_len in sorted(clusters.keys()):
    cluster_IDs = sorted(clusters[CDRH3_len].keys())
    for x in range(len(cluster_IDs)):
      cluster_ID1 = cluster_IDs[x]
      if cluster_ID1 not in clusters[CDRH3_len].keys(): continue
      merge_complete = 0
      while merge_complete == 0:
        merge_count = 0
        for y in range(len(cluster_IDs)):
          cluster_ID2 = cluster_IDs[y]
          if cluster_ID2 not in clusters[CDRH3_len].keys(): continue
          if x < y:
            success_merge = 0
            CDRH3_1s = [CDRH3_dict[Ab1] for Ab1 in clusters[CDRH3_len][cluster_ID1]]
            CDRH3_2s = [CDRH3_dict[Ab2] for Ab2 in clusters[CDRH3_len][cluster_ID2]]
            for CDRH3_1 in CDRH3_1s:
              for CDRH3_2 in CDRH3_2s:
                dist = distance.hamming(CDRH3_1, CDRH3_2)
                if dist/float(len(CDRH3_2)) < cutoff:
                  print ('merging: %s, %s' % (str(CDRH3_len)+'-'+str(cluster_ID1), str(CDRH3_len)+'-'+str(cluster_ID2)))
                  clusters[CDRH3_len][cluster_ID1].extend(clusters[CDRH3_len][cluster_ID2])
                  del clusters[CDRH3_len][cluster_ID2]
                  merge_count =+ 1
                  success_merge = 1
                  break
              if success_merge == 1: break
        if merge_count == 0: merge_complete = 1
  return (clusters)

def writing_cluster_info(clusters, CDRH3_dict, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w') 
  outfile.write("\t".join(map(str,['ID','CDRH3','cluster']))+"\n")
  for CDRH3_len in clusters.keys():
    cluster_count = 0
    for cluster_ID in clusters[CDRH3_len]:
      cluster_count += 1
      for Ab in clusters[CDRH3_len][cluster_ID]:
        outfile.write("\t".join([Ab, CDRH3_dict[Ab], str(CDRH3_len)+'-'+str(cluster_count)])+"\n")
  outfile.close()

def main():
  cutoff   = 0.2
  filename = 'result/CDRH3.tsv'
  outfile  = 'result/CDRH3_cluster.tsv'
  CDRH3_dict = file_to_dict(filename)
  clusters = clustering_CDRH3(CDRH3_dict, cutoff)
  clusters = merging_cluster(clusters, CDRH3_dict, cutoff)
  writing_cluster_info(clusters, CDRH3_dict, outfile)

if __name__ == "__main__":
  main()
