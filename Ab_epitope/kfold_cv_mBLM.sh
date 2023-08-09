#!/bin/bash

# Set the total number of iterations for the loop
total_iterations=15


for ((i=1; i<=$total_iterations; i++)); do
    python script/split_dataset.py
    cp result/epitope_train.tsv dataset/epitope_train_$i.tsv
    cp result/epitope_test.tsv dataset/epitope_test_$i.tsv
    cp result/epitope_val.tsv dataset/epitope_val_$i.tsv
    python script/rm_fas_repeats.py result/epitope_clean.fasta result/epitope_clean_v2.fasta

    python ./Epitope_Clsfr/train.py -n mBLM_attention -lm mBLM -hd 768 -ckn mBLM-$i
    python ./Epitope_Clsfr/test.py -n mBLM_attention -lm mBLM -hd 768 -ckn epoch=10_mBLM-$i.ckpt
    python ./Epitope_Clsfr/predict.py -dp result/Flu_unknown.csv -n mBLM_attention -lm mBLM -hd 768 -ckn epoch=10_mBLM-$i.ckpt 
    python ./Epitope_Clsfr/explain.py -n mBLM_attention -lm mBLM -hd 768 -ckn epoch=10_mBLM-$i.ckpt -o dataset/explain_mBLM_$i/ 
    cp result/mBLM_attention_confusion_matrix.png dataset/mBLM_attention_confusion_matrix$i.png
    cp result/mBLM_attention_epitope_test_prediction.tsv dataset_esm/mBLM_attention_epitope_test_prediction_$i.tsv
    cp result/Flu_unknown_prediction.tsv dataset/Flu_unknown_prediction_$i.tsv
    
    # Calculate the progress percentage
    progress=$((100 * i / total_iterations))

    # Print the progress bar
    echo -ne "Progress: $progress% \r"
    sleep 0.1
done
