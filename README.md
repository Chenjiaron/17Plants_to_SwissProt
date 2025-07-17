# Structure-guided plant protein discovery

## Result dataframe and retrieval system

Results tables and the retrieval system can be found on our website: http://47.111.75.50/ (currently being built).

## Data

Due to the large data volume, the raw data and processed analysis results for this study, as contained within the DATA folder in our code, are also available for download on the website.

## Codes
The code for generating the final results table resides in the `Codes` folder.

The code's operational flow to yield the results table.

1. Stastic basic annotation from EggNog-mapper, Uniprot and Tair:

```
python Codes/get_cluster_all_annotation.py
```

2. Stastic detailed annotation from Tair:

```
python Codes/get_cluster_tair_annotation.py
```

3. Stastic detailed annotation from foldseek alignment result:

```
python Codes/get_cluster_foldseek_annotation.py
```

4. Compare structural alignment results with sequence alignment results and present them in a final table:

```
python Codes/get_cluster_foldseek_blastp_diff.py
```

