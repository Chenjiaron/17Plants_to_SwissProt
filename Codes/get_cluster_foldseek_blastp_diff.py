import pandas as pd
import os
import chardet
from multiprocessing import Pool, cpu_count
import time

# --- Configuration ---
# Determine the number of processes for parallel execution.
# Adjust as needed based on system resources and data volume.
NUM_PROCESSES = min(30, max(1, cpu_count() - 1)) # Use up to 30 processes, or CPU count - 1, whichever is smaller

# --- Functions for Parallel Processing ---
def process_blastp_file(file_path, swiss_entry_ids_set):
    """
    Function to process a single BLASTP result file.
    Args:
        file_path (str): Path to the BLASTP file.
        swiss_entry_ids_set (set): Set of valid UniProt entry IDs for filtering.
    Returns:
        pd.DataFrame or None: Processed DataFrame for the file, or None if an error occurs.
    """
    try:
        # Detect encoding for robust reading.
        with open(file_path, 'rb') as rawdata:
            result = chardet.detect(rawdata.read())
        encoding = result['encoding']

        current_df = pd.read_csv(file_path, header=None, sep='\t', encoding=encoding)

        if current_df.empty:
            return None

        # Ensure enough columns exist.
        if current_df.shape[1] < 12:
            return None

        current_df.columns = [
            'query_id', 'subject_id', 'percent_identity', 'alignment_length',
            'mismatches', 'gap_opens', 'q_start', 'q_end',
            's_start', 's_end', 'e_value', 'bit_score'
        ]

        # Filter by e_value.
        current_df = current_df[current_df['e_value'] < 0.01]

        # Convert to string to prevent errors with .str accessor if column contains non-strings.
        current_df['subject_id'] = current_df['subject_id'].astype(str)
        current_df['query_id'] = current_df['query_id'].astype(str)

        # Extract UniProt Entry IDs using vectorized string operations.
        current_df['subject_id'] = current_df['subject_id'].str.split('|').str[1].fillna('')
        current_df['query_id'] = current_df['query_id'].str.split('|').str[1].fillna('')

        # Filter subject_id based on pre-filtered swiss_entry_ids set.
        current_df = current_df[current_df['subject_id'].isin(swiss_entry_ids_set)]

        # Keep only the best match for each query based on bit_score.
        current_df = current_df.loc[current_df.groupby('query_id')['bit_score'].idxmax()]

        return current_df

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

# --- Main Script Execution ---
start_time = time.time()

# --- 1. Load Initial Data ---
result_df = pd.read_csv('Stastic/Nonsingleton_clusters_allinfo_stastic_add_foldseek.csv')
df_all_swiss_info = pd.read_csv('uniprot_annotate/uniprotkb_reviewed_true_2025_05_29.tsv', sep='\t')

# --- 2. Filter UniProt Entries ---
exclude_terms = [
    'Uncharacterized', 'uncharacterized', 'Uncharacterised', 'uncharacterised',
    'Unknown', 'unknown', 'Hypothetical', 'hypothetical',
    'unannotated', r'DUF\d{2,}', r'putative', r'Putative'
]
exclude_pattern = '|'.join(exclude_terms)

all_info_filtered = df_all_swiss_info[~df_all_swiss_info['Protein names'].str.contains(exclude_pattern, regex=True, na=False)]
swiss_entry_ids = set(all_info_filtered['Entry'])

print(f"Time after initial data loading and UniProt filtering: {time.time() - start_time:.2f} seconds")

# --- 3. Process BLASTP Results in Parallel ---
blastp_results_dir = 'DATA/blastp_result/blastp_reault_all' #Obtained by decompressing DATA/blastp_result_all.tar.gz
blastp_files = [os.path.join(blastp_results_dir, f) for f in os.listdir(blastp_results_dir) if f.endswith('.fasta')]

print(f"Starting parallel processing of {len(blastp_files)} BLASTP files using {NUM_PROCESSES} processes...")
with Pool(NUM_PROCESSES) as pool:
    # Prepare arguments for starmap: (file_path, swiss_entry_ids_set) for each file.
    args = [(f, swiss_entry_ids) for f in blastp_files]
    processed_dfs = pool.starmap(process_blastp_file, args)

# Filter out None values (from files that failed to process or were empty).
processed_dfs = [df for df in processed_dfs if df is not None and not df.empty]

if processed_dfs:
    combined_blastp_df = pd.concat(processed_dfs, ignore_index=True)
else:
    print("No valid BLASTP results were processed successfully.")
    # Create an empty DataFrame with expected columns if no files were processed.
    combined_blastp_df = pd.DataFrame(columns=[
        'query_id', 'subject_id', 'percent_identity', 'alignment_length',
        'mismatches', 'gap_opens', 'q_start', 'q_end',
        's_start', 's_end', 'e_value', 'bit_score'
    ])

print(f"Time after parallel BLASTP file processing: {time.time() - start_time:.2f} seconds")

# --- 4. Build blastp_results_dict Efficiently ---
if not combined_blastp_df.empty:
    # Group by query_id and aggregate subject_ids into lists.
    blastp_results_dict = combined_blastp_df.groupby('query_id')['subject_id'].apply(list).to_dict()
else:
    blastp_results_dict = {}
    print("Combined BLASTP DataFrame is empty, blastp_results_dict will be empty.")

print(f"Time after building blastp_results_dict: {time.time() - start_time:.2f} seconds")

# --- 5. Determine 'blastp_check' Column Using Vectorized Methods ---
# Ensure required columns exist in result_df.
if 'target' not in result_df.columns:
    raise ValueError("Error: 'target' column not found in result_df. Please check your data and adjust column name.")
if 'represent' not in result_df.columns:
    raise ValueError("Error: 'represent' column not found in result_df. Please check your data and adjust column name.")

# Convert 'represent' and 'target' columns to string type for consistent operations.
result_df['represent'] = result_df['represent'].astype(str)
result_df['target'] = result_df['target'].astype(str)

# Map 'represent' column to its corresponding list of subject_ids from the dictionary.
# Use .get(key, []) to handle cases where 'represent' ID might not be in blastp_results_dict.
result_df['blastp_hits_for_represent'] = result_df['represent'].apply(
    lambda x: blastp_results_dict.get(x, [])
)

# Check if 'target' is present in the mapped list of hits.
result_df['blastp_check'] = result_df.apply(
    lambda row: '1' if row['target'] in row['blastp_hits_for_represent'] else '0',
    axis=1
)

# Drop the temporary column as it's no longer needed.
result_df = result_df.drop(columns=['blastp_hits_for_represent'])

print(f"Total execution time: {time.time() - start_time:.2f} seconds")

print("\n'blastp_check' column successfully added to result_df:")
print(result_df[['represent', 'target', 'blastp_check']].head())

# Save the modified result_df.
result_df.to_csv('Stastic/Nonsingleton_clusters_allinfo_stastic_add_foldseek_with_blastp_check.csv', index=False)