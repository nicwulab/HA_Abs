#!/bin/bash

# Set the total number of iterations for the loop
total_iterations=40


for ((i=1; i<=$total_iterations; i++)); do
    python script/split_dataset.py
    cp result/epitope_train.tsv dataset_esm/epitope_train_$i.tsv
    cp result/epitope_test.tsv dataset_esm/epitope_test_$i.tsv
    cp result/epitope_val.tsv dataset_esm/epitope_val_$i.tsv
    python script/rm_fas_repeats.py result/epitope_clean.fasta result/epitope_clean_v2.fasta

    python ./Epitope_Clsfr/train.py -n esm2_attention -lm esm2_t33_650M_UR50D -hd 1280 -ckn esm2-$i
    python ./Epitope_Clsfr/test.py -n esm2_attention -lm esm2_t33_650M_UR50D -hd 1280 -ckn epoch=29_esm2-$i.ckpt
    python ./Epitope_Clsfr/predict.py -dp result/Flu_unknown.csv-n esm2_attention -lm esm2_t33_650M_UR50D -hd 1280 -ckn epoch=29_esm2-$i.ckpt 
    python ./Epitope_Clsfr/explain.py -n esm2_attention -lm esm2_t33_650M_UR50D -hd 1280 -ckn epoch=29_esm2-$i.ckpt -o dataset_esm/explain_esm2_$i/ 
    cp result/esm2_attention_confusion_matrix.png dataset_esm/esm2_attention_confusion_matrix_$i.png
    cp result/esm2_attention_epitope_test_prediction.tsv dataset_esm/esm2_attention_epitope_test_prediction_$i.tsv
    cp result/Flu_unknown_prediction.tsv dataset_esm/Flu_unknown_prediction_$i.tsv
    
    # Calculate the progress percentage
    progress=$((100 * i / total_iterations))

    # Print the progress bar
    echo -ne "Progress: $progress% \r"
    sleep 0.1
done
