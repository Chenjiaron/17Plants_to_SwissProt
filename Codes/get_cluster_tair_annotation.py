import pandas as pd
import re
import os

# --- Load DataFrames ---
result_df = pd.read_csv('Stastic/Nonsingleton_clusters_allinfo_stastic.csv')
tair_pep_file = 'DATA/uni2tair/Araport11_pep_20250214'
tair_description_file = 'DATA/tair_annotate/functional_description/Araport11_functional_descriptions_20240331.csv'
uni2tair_file = 'DATA/uni2tair/Uniprot2AGI.txt'
represent2tair_file = 'DATA/represent2tair/results.csv'
df_all_swiss_info = pd.read_csv('DATA/uniprot_annotate/uniprotkb_reviewed.tsv', sep='\t')

# Select and rename columns from UniProt info.
df_swiss_info = df_all_swiss_info[['Entry', 'Protein names', 'Organism']]
df_swiss_info.columns = ['Entry', 'Protein names(uniprotkb)(seq)', 'target_specie(seq)']

# Load and process represent2tair_file (BLAST results for UniProt to TAIR).
re2tair_df = pd.read_csv(represent2tair_file, header=None, sep='\t')
re2tair_df.columns = [
    'query_id', 'subject_id', 'percent_identity', 'alignment_length',
    'mismatches', 'gap_opens', 'q_start', 'q_end',
    's_start', 's_end', 'e_value', 'bit_score'
]
re2tair_df = re2tair_df[re2tair_df['e_value'] < 0.01]
re2tair_df['query_id'] = re2tair_df['query_id'].str.split('|').str[1] # Extract UniProt ID.

# --- Define exclusion patterns for annotations ---
exclude_terms = [
    'Uncharacterized', 'uncharacterized', 'Uncharacterised', 'uncharacterised',
    'Unknown', 'unknown', 'Hypothetical', 'hypothetical',
    'unannotated', r'DUF\d{2,}', r'putative', r'Putative'
]
exclude_pattern = '|'.join(exclude_terms)

# Filter UniProt data based on exclusion pattern.
all_info = df_all_swiss_info[~df_swiss_info['Protein names(uniprotkb)(seq)'].str.contains(exclude_pattern, regex=True, na=False)]
swiss_entry_ids = set(all_info['Entry'])

# Load TAIR functional descriptions.
df_tair_des = pd.read_csv(tair_description_file)
df_tair_des = df_tair_des[['name', 'Curator_summary', 'Computational_description']]
des_dict_cur = df_tair_des.set_index('name')['Curator_summary'].to_dict()
des_dict_com = df_tair_des.set_index('name')['Computational_description'].to_dict()

# --- Function to check if an annotation is valid ---
def is_valid_annotation(annotation):
    if pd.isna(annotation) or str(annotation).strip() == "":
        return False
    if re.search(exclude_pattern, str(annotation), re.IGNORECASE):
        return False
    return True

# --- Process BLASTP results from directory ---
blastp_results_dir = 'blastp_result'
all_blastp_dfs = []

for filename in os.listdir(blastp_results_dir):
    if filename.endswith('.fasta'): # Assuming BLASTP output files are .fasta (or adjust extension).
        file_path = os.path.join(blastp_results_dir, filename)
        current_df = pd.read_csv(file_path, header=None, sep='\t')
        current_df.columns = [
            'query_id', 'subject_id', 'percent_identity', 'alignment_length',
            'mismatches', 'gap_opens', 'q_start', 'q_end',
            's_start', 's_end', 'e_value', 'bit_score'
        ]
        current_df = current_df[current_df['e_value'] < 0.01]
        current_df['target_species'] = current_df['subject_id'].str.split('_').str[1]
        current_df['subject_id'] = current_df['subject_id'].str.split('|').str[1] # Extract UniProt ID from subject.
        current_df['query_id'] = current_df['query_id'].str.split('|').str[1] # Extract UniProt ID from query.
        current_df = current_df[current_df['subject_id'].isin(swiss_entry_ids)] # Filter by valid UniProt entries.
        # Keep only the best match for each query based on bit_score.
        current_df = current_df.loc[current_df.groupby('query_id')['bit_score'].idxmax()]
        all_blastp_dfs.append(current_df)

# Concatenate all BLASTP DataFrames.
blastp_result_df = pd.concat(all_blastp_dfs, ignore_index=True)

# --- Create AT to Symbol mapping from TAIR pep file ---
at_to_symbol = {}
with open(tair_pep_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            match = re.search(r'>(AT\dG\d{5}\.\d) \| Symbols: (.+?) \|', line)
            if match:
                at_id = match.group(1)
                symbols = match.group(2)
                symbol_list = [symbol.strip() for symbol in symbols.split(',')] if symbols else []
                at_to_symbol[at_id] = symbol_list

# --- Create UniProt to AGI mapping ---
uni_to_agi = {}
with open(uni2tair_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if ';' in line: # Handle cases where AGI IDs are semicolon-separated in the file
            uni_id, agi_id = parts
            uni_to_agi[uni_id] = [agi_id]
        else:
            uni_id, agi_ids_str = parts[0], parts[1]
            agi_ids_list = agi_ids_str.split(';')
            uni_to_agi[uni_id] = agi_ids_list

# --- Add TAIR and AGI information to result_df ---
result_df_extended = result_df.copy()

for index, row in result_df.iterrows():
    uni_id = row['represent']
    agi_ids = uni_to_agi.get(uni_id, [])

    if agi_ids:
        description_curator = []
        description_computational = []
        symbol_ids = []
        for agi_id in agi_ids:
            if agi_id in at_to_symbol:
                symbol_ids.extend(at_to_symbol[agi_id])
            if agi_id in des_dict_cur and is_valid_annotation(des_dict_cur[agi_id]):
                description_curator.append(str(des_dict_cur[agi_id]))
            if agi_id in des_dict_com and is_valid_annotation(des_dict_com[agi_id]):
                description_computational.append(str(des_dict_com[agi_id]))

        result_df_extended.at[index, 'tair_symbol'] = '; '.join(symbol_ids) if symbol_ids else None
        result_df_extended.at[index, 'tair_description'] = '; '.join(description_curator) if description_curator else None
        result_df_extended.at[index, 'tair_description(computational)'] = '; '.join(description_computational) if description_computational else None
        result_df_extended.at[index, 'AGI_id'] = '; '.join(agi_ids) if agi_ids else None

# Drop the second column (index 1) from the extended DataFrame.
result_df_final = result_df_extended.drop(result_df_extended.columns[1], axis=1)

# --- Merge with BLASTP results and UniProt info ---
merged_df = pd.merge(result_df_final, blastp_result_df, left_on='represent', right_on='query_id', how='left')
final_df = merged_df.drop(columns=['query_id'])
final_df_with_swiss_info = pd.merge(final_df, df_swiss_info, left_on='subject_id', right_on='Entry', how='left')
final_df_with_swiss_info = final_df_with_swiss_info.drop(columns=['Entry'])

# Display and save the final DataFrame.
print(final_df_with_swiss_info)
final_df_with_swiss_info.to_csv('Stastic/Nonsingleton_clusters_allinfo_stastic_add_tair_blastp.csv', index=False)