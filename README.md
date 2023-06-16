# Sequence analysis of influenza hemagglutinin (HA) antibodies

## Input files
* [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx): List of HA antibodies collected from publications and patents
* [./doc/all_paired_antibodies_from_GB_v6.xlsx](./doc/all_paired_antibodies_from_GB_v6.xlsx): List of HA antibodies collected from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)

## CDR H3 analysis
1. Extract CDR H3 sequences and references   
``python3 script/parse_Ab_table.py``
    - Input file:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
    - Output files:
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
      - [./result/refs.txt](./result/refs.txt)

2. Clustering CDR H3 sequences   
``python3 script/CDRH3_clustering_optimal.py``
    - Input file:
      - [./result/CDRH3.tsv](./result/CDRH3.tsv)
    - Output file:
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)

3. Analyzing CDR H3 clustering results   
``python3 script/analyze_CDRH3_cluster.py``
    - Input files:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)
    - Output files:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
      - [./result/CDRH3_cluster_summary.tsv](./result/CDRH3_cluster_summary.tsv)

4. Analyzing CDR H3 property   
``python3 script/analyze_CDRH3_property.py``
    - Input files:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
      - [./doc/all_paired_antibodies_from_GB_v6.xlsx](./doc/all_paired_antibodies_from_GB_v6.xlsx)
      - [./result/HA_Abs_clonotype.xlsx](./result/HA_Abs_clonotype.xlsx)
    - Ouput files:
      - [./result/CDRH3_property.tsv](./result/CDRH3_property.tsv)
      - [./result/Ab_for_model.tsv](./result/Ab_for_model.tsv)

5. Create sequence logos for different CDR H3 clusters   
``python3 script/CDRH3_seqlogo.py``
    - Input file:
      - [./result/Ab_info_CDRH3_clustering.tsv](./result/Ab_info_CDRH3_clustering.tsv)
    - Output file:
      - ./CDRH3_seqlogo/*.png

6. Plot CDR H3 property for HA head and stem antibodies   
``Rscript script/plot_CDRH3_property.R``
    - Input file:
      - [./result/CDRH3_property.tsv](./result/CDRH3_property.tsv)
    - Output files:
      - [./graph/CDRH3_length.png](./graph/CDRH3_length.png)
      - [./graph/CDRH3_hydrophobicity.png](./graph/CDRH3_hydrophobicity.png)
      - [./graph/tip_hydrophobicity.png](./graph/tip_hydrophobicity.png)

## Germline usage analysis
1. Clonotype assignment   
``python3 script/assign_clonotype.py``
    - Input files:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
      - [./result/CDRH3_cluster.tsv](./result/CDRH3_cluster.tsv)
    - Output file:
      - [./result/HA_Abs_clonotype.xlsx](./result/HA_Abs_clonotype.xlsx)

2. Compute germline usag and extract public clonotype   
``python3 script/extract_public_clonotype_VDJ.py``
    - Input files:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
      - [./doc/all_paired_antibodies_from_GB_v6.xlsx](./doc/all_paired_antibodies_from_GB_v6.xlsx)
      - [./result/HA_Abs_clonotype.xlsx](./result/HA_Abs_clonotype.xlsx)
    - Output files:
      - [./result/IGHV_freq.tsv](./result/IGHV_freq.tsv)
      - [./result/IGHD_freq.tsv](./result/IGHD_freq.tsv)
      - [./result/IGLV_freq.tsv](./result/IGLV_freq.tsv)
      - [./result/IGV_pair_freq.tsv](./result/IGV_pair_freq.tsv)
      - [./result/HA_Abs_clonotype_public.xlsx](./result/HA_Abs_clonotype_public.xlsx)
      
3. Extract IGHD4-17-encoded head antibodies 
``python3 script/analyze_IGHD4-17.py``   
    - Input file:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
    - Output file:
      - [./result/Abs_IGHD4-17.tsv](./result/Abs_IGHD4-17.tsv)

4. Analyzing the occurrence of YGD motif in CDR H3   
``python3 script/analyze_YGD_motif.py``
    - Input files:
      - [./doc/HA_Abs_v17.xlsx](./doc/HA_Abs_v17.xlsx)
      - [./doc/all_paired_antibodies_from_GB_v6.xlsx](./doc/all_paired_antibodies_from_GB_v6.xlsx)
    - Ouput file:
      - [./result/YGD_motif_freq.tsv](./result/YGD_motif_freq.tsv)

5. Plot VDJ gene usage   
``Rscript script/plot_VDJgene_freq.R``
    - Input files:
      - [./result/IGHV_freq.tsv](./result/IGHV_freq.tsv)
      - [./result/IGHD_freq.tsv](./result/IGHD_freq.tsv)
      - [./result/IGLV_freq.tsv](./result/IGLV_freq.tsv)
    - Output files:
      - [./graph/IGHV_usage.png](./graph/IGHV_usage.png)
      - [./graph/IGHD_usage.png](./graph/IGHD_usage.png)
      - [./graph/IGLV_usage.png](./graph/IGLV_usage.png)

6. Plot IGHV/IGK(L)V pairing frequency   
``Rscript script/plot_Vpair_heatmap.R``
    - Input files:
      - [./result/IGV_pair_freq.tsv](./result/IGV_pair_freq.tsv)
    - Output file:
      - [./graph/Vpair_heatmap.png](./graph/Vpair_heatmap.png)

7. Plot frequency of YGD motif   
``Rscript script/plot_YGD_freq.R``
    - Input file:
      - [./result/YGD_motif_freq.tsv](./result/YGD_motif_freq.tsv)
    - Output file:
      - [./graph/YGD_motif_freq.png](./graph/YGD_motif_freq.png)
