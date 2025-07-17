import pandas as pd

# --- Load DataFrames ---
# Assuming these files and paths are correct and accessible.
result_df = pd.read_csv('Stastic/Nonsingleton_clusters_allinfo_stastic_add_tair_blastp.csv')
df_foldseek_result = pd.read_csv('DATA/foldseek_cluster_result/plant2swiss_v2_bestmatches_all.csv')
df_foldseek_result['query'] = df_foldseek_result['query'].str.split('-').str[1]
df_foldseek_result['target'] = df_foldseek_result['target'].str.split('-').str[1]

# Load UniProt information.
df_all_swiss_info = pd.read_csv('DATA/uniprot_annotate/uniprotkb_reviewed.tsv', sep='\t')
df_tmalign = pd.read_csv('DATA/TM_align_result/tm_result.txt',sep = '\t')

# Assign columns to tmalign DataFrame.
df_tmalign.columns = ['represent','target','tmscore_tmalign','seq_id','rmsd','aligned_length']

# Prepare result_df for merging.
result_df_cleaned = result_df.copy()
result_df_cleaned = result_df_cleaned.drop(columns=['target_species'])

# Select and rename columns from UniProt info for merging.
df_swiss_info_subset = df_all_swiss_info[['Entry','Protein names','Organism']]
df_swiss_info_subset.columns = ['Entry','Protein names(uniprotkb)(str)','target_specie(str)']

df_swiss_info_protein_names = df_all_swiss_info[['Entry','Protein names']]
df_swiss_info_protein_names.columns = ['Entry','Protein names(uniprotkb)(reviewed)']

# --- Merge DataFrames ---
# Merge with foldseek results.
merged_df = pd.merge(result_df_cleaned, df_foldseek_result, left_on='represent', right_on='query', how='left')
final_df = merged_df.drop(columns=['query'])

# Merge with UniProt protein names (reviewed).
final_df_uniprot1 = pd.merge(final_df, df_swiss_info_protein_names, left_on='reviewed_protein', right_on='Entry', how='left')
final_df_uniprot1 = final_df_uniprot1.drop(columns=['Entry'])

# Merge with UniProt protein names and organism.
final_df_uniprot2 = pd.merge(final_df_uniprot1, df_swiss_info_subset, left_on='target', right_on='Entry', how='left')
final_df_uniprot2 = final_df_uniprot2.drop(columns=['Entry'])

# Merge with TM-align results.
final_df_combined = pd.merge(final_df_uniprot2, df_tmalign, on=['represent','target'], how='left')

# Save the final DataFrame to CSV.
final_df_combined.to_csv('Stastic/Nonsingleton_clusters_allinfo_stastic_add_foldseek.csv', index=False)