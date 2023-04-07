This README describes fine-tune ESM protein language model for antibody binding specificity prediction

# Dataset clean and split
Firstly, we remove "unknown" group sequence that is similar to annotated groups (sequence identity 0.9), and then define "unknown" as "Others"

Then, we split into train/test/val by sequence identity at maximum 80%
```
python split_dataset.py
```
# ESM embbeding extraction
```
python rm_fas_repeats.py result/epitope_clean.fasta result/epitope_clean_v2.fasta
python esm_extractor.py esm2_t33_650M_UR50D result/epitope_clean_v2.fasta result/esm2_t33_650M_embedding --repr_layers 33 --include per_tok contacts
```
# Model build
### fine-tune language model for classification

[model.py](./epitope_classifier/model.py)


# Model training
### Benchmarking
```
python train.py -n onehot_baseline -lm onehot -hd 1280
python train.py -n esm2_baseline -lm esm2_t33_650M_UR50D -hd 1280
```

# Model testing
### test model performance on test set and output a 'epitope_test_prediction.tsv' and confusion matrix

```
python test.py -n onehot_baseline -lm onehot -hd 1280 -ckn epoch=10-step=4202.ckpt
python test.py -n esm2_baseline -lm esm2_t33_650M_UR50D -hd 1280 -ckn epoch=10-step=4202.ckpt
```
# Model predicting for new dataset
### predict on new dataset with input parameter '-dp' for dataframe path
```
python predict.py -dp /home/yiquan2/ESM_Ab/HA_Abs/HA_epitope/result/Flu_unknown.csv -n esm2_baseline -lm esm2_t33_650M_UR50D -hd 1280 -ckn epoch=10-step=4202.ckpt
```