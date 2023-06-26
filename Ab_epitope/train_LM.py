from transformers import RobertaConfig, RobertaForMaskedLM
from transformers import AutoTokenizer,Trainer,  TrainingArguments, DataCollatorForLanguageModeling
from torch.utils.data import Dataset
from transformers import PreTrainedTokenizer
from typing import Dict
import torch
import math
import os
import pandas as pd
import argparse

class LineByLineTextDataset(Dataset):


    def __init__(self, tokenizer: PreTrainedTokenizer, file_path: str, block_size: int):
        
        if os.path.isfile(file_path) is False:
            raise ValueError(f"Input file path {file_path} not found")

        with open(file_path, encoding="utf-8") as f:
            lines = [line for line in f.read().splitlines() if (len(line) > 0 and not line.isspace())]

        batch_encoding = tokenizer(lines, add_special_tokens=True, truncation=True, max_length=block_size)
        self.examples = batch_encoding["input_ids"]
        self.examples = [{"input_ids": torch.tensor(e, dtype=torch.long)} for e in self.examples]

    def __len__(self):
        return len(self.examples)

    def __getitem__(self, i) -> Dict[str, torch.tensor]:
        return self.examples[i]
    
def prepare_dataset(train_data,val_data,lm_data_dir):
    # load dataset
    if train_data.split('.')[-1] == 'tsv':
        train_df = pd.read_csv(train_data,sep='\t')
        val_df = pd.read_csv(val_data,sep='\t')
    elif train_data.split('.')[-1] == 'xlsx':
        train_df = pd.read_excel(train_data)
        val_df = pd.read_excel(val_data)
    elif train_data.split('.')[-1] == 'csv':
        train_df = pd.read_csv(train_data)
        val_df = pd.read_csv(val_data)
    else:
        print('ERROR: unfined data format')
        
    train_sequences = train_df['sequence_alignment_aa_heavy'].apply(lambda x: ' '.join(list(x))) + ' <null_1> ' + train_df['sequence_alignment_aa_light'].apply(lambda x: ' '.join(list(x)))
    
    val_sequences = val_df['sequence_alignment_aa_heavy'].apply(lambda x: ' '.join(list(x))) + ' <null_1> ' + val_df['sequence_alignment_aa_light'].apply(lambda x: ' '.join(list(x)))

    with open(os.path.join(lm_data_dir,'LM_train.txt') , 'w') as f:
        for item in train_sequences.tolist():
            f.write("%s\n" % item)
    with open(os.path.join(lm_data_dir,'LM_val.txt') , 'w') as f:
        for item in val_sequences.tolist():
            f.write("%s\n" % item)
            
    return os.path.join(lm_data_dir,'LM_train.txt'), os.path.join(lm_data_dir,'LM_val.txt')

def mBLM(vocab_size,max_position_embeddings,num_attention_heads,num_hidden_layers):
    
    config = RobertaConfig(
        vocab_size=vocab_size,
        max_position_embeddings=max_position_embeddings,
        num_attention_heads=num_attention_heads,
        num_hidden_layers=num_hidden_layers,
        type_vocab_size=1,
    )


    model = RobertaForMaskedLM(config=config)
    
    return model

def train_mBLM(tokenizer,train_dir,val_dir,epoch,per_gpu_train_batch_size,output_dir):
    train_dataset = LineByLineTextDataset(
        tokenizer=tokenizer,
        file_path=train_dir,
        block_size=150,
    )
    val_dataset = LineByLineTextDataset(
        tokenizer=tokenizer,
        file_path=val_dir,
        block_size=150,
    )

    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer, mlm=True, mlm_probability=0.15
    )

    training_args = TrainingArguments(
        output_dir=output_dir,
        overwrite_output_dir=True,
        num_train_epochs=epoch,
        per_gpu_train_batch_size=per_gpu_train_batch_size,
        save_steps=1000,
        save_total_limit=2,
        prediction_loss_only=True,
    )

    trainer = Trainer(
        model = model,
        args = training_args,
        data_collator = data_collator,
        train_dataset = train_dataset,
        eval_dataset = val_dataset
    )

    trainer.train()
    eval_results = trainer.evaluate()
    print(f"Evaluation perplexity: {math.exp(eval_results['eval_loss']):.2f}")
    # save model
    trainer.save_model(args.model_location)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train custom antibody language model")
    parser.add_argument(
        "--model_location",
        default='./mBLM',
        type=str,
        help="model dir to save ",)
    parser.add_argument(
        "--tokenizer",
        default='facebook/esm2_t6_8M_UR50D',
        type=str,
        help="tokenizer name to use",)
    parser.add_argument(
        "--vocab_size",
        default=33,
        type=int,
        help=" the token size used in the language model ",)
    parser.add_argument(
        "--max_position_embeddings",
        default=514,
        type=int,
        help=" the maximum sequence length that the model can accept as input ",)
    parser.add_argument(
        "--num_attention_heads",
        default=12,
        type=int,
        help="the number of attention head ",)
    parser.add_argument(
        "--num_hidden_layers",
        default=6,
        type=int,
        help="the number of hidden layer ",)
    parser.add_argument(
        "--epoch",
        default=20,
        type=int,
        help="the number of training epoch ",)
    parser.add_argument(
        "--train_data",
        default='./data/dataset/OAS_memory_val.tsv',
        type=str,
        help="train dataset",)
    parser.add_argument(
        "--val_data",
        default='./data/dataset/OAS_memory_val.tsv',
        type=str,
        help="validation dataset",)
    parser.add_argument(
        "--lm_data_dir",
        default='./data/dataset/',
        type=str,
        help="Language model dataset directory ",)
    parser.add_argument(
        "--batch_size", 
        type=int, default=64, 
        help="batch size per GPU")
    args = parser.parse_args()
    
    # load tokenizer
    os.makedirs(args.model_location, exist_ok=True)
    tokenizer = AutoTokenizer.from_pretrained(args.tokenizer)
    tokenizer.save_pretrained(args.model_location)

    # load dataset
    train_dir, val_dir = prepare_dataset(args.train_data,args.val_data,args.lm_data_dir)
    
    # define model
    model = mBLM(vocab_size=args.vocab_size,
                 max_position_embeddings=args.max_position_embeddings,
                 num_attention_heads=args.num_attention_heads,
                 num_hidden_layers=args.num_hidden_layers)
    # train model
    train_mBLM(tokenizer=tokenizer,
               train_dir=train_dir,
               val_dir=val_dir,
               epoch=args.epoch,
               per_gpu_train_batch_size=args.batch_size,
               output_dir=args.model_location)
    

