# 17Plants_to_SwissProt

## Result dataframe and retrieval system

Results tables and the retrieval system can be found on our website: http://47.111.75.50/ (currently being built).

## Data

The raw data and process analysis results for this study are all contained within the `DATA` folder.

## Codes
The code for generating the final results table resides in the `Codes` folder.

The code's operational flow to yield the results table.

1. Clone the repository:

```
git clone https://github.com/westlake-repl/ESM-Ezy.git
```

2. Install the required packages:

```
conda env create -f environment.yml
```

3. Download the pre-trained ESM-1b model:

```
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt -O ckpt/esm1b_t33_650M_UR50S.pt
wget https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt -O ckpt/esm1b_t33_650M_UR50S-contact-regression.pt
```

4. Train ESM-Ezy:

```
python scripts/train.py --train_positive_data data/train/train_positive.fa --train_negative_data data/train/train_negative.fa --test_positive_data data/train/test_positive.fa --test_negative_data data/train/test_negative.fa --model_path ckpt/esm1b_t33_650M_UR50S.pt
```

