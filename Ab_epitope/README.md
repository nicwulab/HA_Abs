# Memory-B-cells Language Model (mBLM) for antibody epitope prediction guidelines
This README describes the mBLM in the paper: [paper link]()


## Contents

* [Env setup](#Env-setup)   
* [Dataset download and process](#dataset-download-and-process) 
* [Train memory B cell Language Model](#train-memory-b-cell-language-model)
* [Epitope prediction](#epitope-prediction)   
* [Antibody binding sites identification](#antibody-binding-sites-identification)   


## Env setup

### if you set up env using conda, run conda installation as follow:
```commandline
conda env create -f environment.yml
```

## Dataset download and pre-process

### Dataset

- OSA human paired memory B cell: [https://opig.stats.ox.ac.uk/webapps/oas/oas_paired/](https://opig.stats.ox.ac.uk/webapps/oas/oas_paired/)
- Flu Antibody dataset in this paper: [Flu](raw_data/)
- SARS-CoV-2 Antibody dataset in this paper: [SARS-CoV-2](raw_data/)

### dataset for memory-B-cells Language Model
We downloaded and processed all OAS memory paired B cell seuqences from OAS.
Then, heavy chain was clustered by 95% sequence identity using [cdhit](https://github.com/weizhongli/cdhit). Sequence representative was chosen from each cluster.
For training and test purpose, the dataset was splitted by different sequence identity (50%, 60%, 70%, 80%, 90%).

```commandline
./data_clean_for_LM.sh
```
### dataset for antibody epitope prediction
Note: for epitope prediction, heavy chain sequence was used only. 
Firstly, we remove "unknown" sequence that is similar to annotated sequence (sequence identity 0.9), and then define "unknown" as "Others".
Then, we split into train/test/val by sequence identity at maximum 80% (>= 26 AA differences).

```commandline
./data_clean_for_epitope.sh
```

## Train memory B cell Language Model

mBLM was adapted from RoBERTa model [RoBERTa: A Robustly Optimized BERT Pretraining Approach](https://arxiv.org/abs/1907.11692).
model and training details see paper[paper link]().

```commandline
python train_LM.py
```

## Epitope prediction
we then fine-tuned mBLM for multi-epitopes prediction.

### extract mBLM embedding features

```commandline
python extract_mBLM_feature.py --model_location mBLM --fasta_file result/epitope_info_clstr_v2.fasta --output_dir result/mBLM_embedding
```
### benchmarking
```commandline
python ./Epitope_Clsfr/train.py -n onehot_baseline -lm onehot -hd 26
python ./Epitope_Clsfr/train.py -n esm2_attention -lm esm2_t33_650M_UR50D -hd 1280
python ./Epitope_Clsfr/train.py -n mBLM_attention -lm mBLM -hd 768
```
### model test and predict
```commandline
python ./Epitope_Clsfr/test.py -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM.ckpt

python ./Epitope_Clsfr/predict.py -dp result/Flu_unknown.csv -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM.ckpt

python ./Epitope_Clsfr/predict.py -dp result/Sarah_stem_antibodies.xlsx -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM.ckpt
```

## Antibody binding sites identification
In order to investigate how model make a accurate prediction, Grad-CAM [https://arxiv.org/abs/1610.02391](https://arxiv.org/abs/1610.02391) technique was used to show model capturing binding sites in a residues level.

### Grad-CAM
```commandline
python ./Epitope_Clsfr/explain.py -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM.ckpt

python ./Epitope_Clsfr/explain.py -dfp result/Flu_unknown.csv -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM.ckpt -o result/Flu_unknown_explain/ --provide_dataset
```

### Visualization
python script/write_pdb_visualizer.py

python script/pymol_viz.py

python script/cluster_saliency_map.py