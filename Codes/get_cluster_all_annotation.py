import pandas as pd
import re, os

# File paths
blastp_file = 'DATA/blastp_result/blast_results.csv'
eggnogmapper_result = 'DATA/eggnog_result/final_cluster.emapper.annotations.tsv'
cluster_file = 'DATA/foldseek_cluster_result/foldseek_cluster_result.txt'
uni2tair_file = 'DATA/uni2tair/uni2tair_blastp_result.txt'
tair_annotate_dir = 'DATA/tair_annotate/'

# Load pLDDT scores
with open('uniprot_annotate/Plant_structure_all_plddt.txt', 'r') as file:
    plddt_dict = {parts[0]: float(parts[1]) for line in file if (parts := line.strip().split()) and len(parts) == 2}


def find_highest_value_item(item_list, value_dict):
    highest_item = None
    highest_value = float('-inf')

    for item in item_list:
        if item in value_dict:
            current_value = value_dict[item]
            if isinstance(current_value, (int, float)) and current_value > highest_value:
                highest_value = current_value
                highest_item = item
    return highest_item


# --- Process BLASTP results ---
# This block reads multiple BLASTP result files, filters them by e-value,
# extracts UniProt IDs, removes duplicates, and concatenates them into a single DataFrame.
blastp_results_dir = 'blastp_result'
all_blastp_dfs = []
for filename in os.listdir(blastp_results_dir):
    if filename.endswith('.fasta'):
        file_path = os.path.join(blastp_results_dir, filename)
        current_df = pd.read_csv(file_path, header=None, sep='\t')
        current_df.columns = [
            'query_id', 'subject_id', 'percent_identity', 'alignment_length',
            'mismatches', 'gap_opens', 'q_start', 'q_end',
            's_start', 's_end', 'e_value', 'bit_score'
        ]
        current_df = current_df[current_df['e_value'] < 0.01]
        current_df['query_id'] = current_df['query_id'].str.split('|').str[1]
        current_df.drop_duplicates(subset=['query_id'], inplace=True)
        all_blastp_dfs.append(current_df)
blastp_result_df = pd.concat(all_blastp_dfs, ignore_index=True)

# --- Process EggNOG-mapper results ---
# This section loads EggNOG annotations, filters by e-value, extracts UniProt IDs,
# and pre-compiles a regex pattern for excluding uninformative annotations.
eggnog_result_df = pd.read_csv(eggnogmapper_result, header=4, sep='\t')
eggnog_result_df = eggnog_result_df[eggnog_result_df['evalue'] < 0.01]
eggnog_result_df['#query'] = eggnog_result_df['#query'].str.split('|').str[1]
exclude_pattern = re.compile(
    '|'.join([
        'Uncharacterized', 'uncharacterized', 'Uncharacterised', 'uncharacterised',
        'Unknown', 'unknown', 'Hypothetical', 'hypothetical',
        'unannotated', r'DUF\d{2,}', r'putative', r'Putative'
    ]), re.IGNORECASE
)

# --- Load Foldseek clusters ---
# Reads the cluster file and filters for clusters containing more than one element, then assigns names.
clusters = []
with open(cluster_file, 'r') as file:
    for line in file:
        elements = line.strip().split(' ')
        if len(elements) > 1:
            clusters.append(elements)
named_clusters = {f'cluster_{i + 1}': cluster for i, cluster in enumerate(clusters)}

# --- Load UniProt annotation data ---
# This block loads additional UniProt information (reviewed status and species name)
# and creates dictionaries for quick lookups.
df_all_info_review = pd.read_csv('uniprot_annotate/uniprotkb_reviewed_true.tsv', sep='\t')
df_all_info = pd.read_excel('uniprot_annotate/all_plant_structures_info_from_uniprot.xlsx')

entry_info_spe = {row['Entry']: {'spe_name': row['spe_name']} for _, row in df_all_info.iterrows()}
entry_info_reviewed = {row['Entry']: {'Reviewed': row['Reviewed']} for _, row in df_all_info_review.iterrows()}
entry_info_egg = {
    row['#query']: {'description': row['Description'], 'pfam': row['PFAMs'], 'preferred_name': row['Preferred_name']}
    for _, row in eggnog_result_df.iterrows()}
entry_info_blastp = {row['query_id']: {'target': row['subject_id']} for _, row in blastp_result_df.iterrows()}

# --- Summarize cluster information ---
# Iterates through each named cluster, extracts various annotations (species, reviewed status,
# EggNOG annotations, BLASTP results) for its constituent proteins, and aggregates counts.
cluster_summary = {}
for cluster_name, elements in named_clusters.items():
    summary_data = {
        'spe_name': {},
        'reviewed_count': 0,
        'reviewed_name': [],
        'total_count': len(elements),
        'eggnog_annotation_count': 0,
        'eggnog_pfam_count': 0,
        'eggnog_pname_count': 0,
        'blastp_result_count': 0
    }

    for element in elements:
        entry_id = element.split('-')[1]
        reviewed = entry_info_reviewed.get(entry_id, {}).get('Reviewed', '')
        spe_info = entry_info_spe.get(entry_id, {})
        spe_name = spe_info.get('spe_name', '')

        egg_info = entry_info_egg.get(entry_id, {})
        egg_description = str(egg_info.get('description', ''))
        egg_pfam = str(egg_info.get('pfam', ''))
        egg_pname = str(egg_info.get('preferred_name', ''))

        blastp_result = entry_info_blastp.get(entry_id, {}).get('target', '')

        if spe_name:
            summary_data['spe_name'].setdefault(spe_name, []).append(element)

        if reviewed == 'reviewed':
            summary_data['reviewed_count'] += 1
            summary_data['reviewed_name'].append(entry_id)

        if not exclude_pattern.search(egg_description) and egg_description != '-':
            summary_data['eggnog_annotation_count'] += 1
        if not exclude_pattern.search(egg_pfam) and egg_pfam != '-':
            summary_data['eggnog_pfam_count'] += 1
        if not exclude_pattern.search(egg_pname) and egg_pname != '-':
            summary_data['eggnog_pname_count'] += 1
        if blastp_result:
            summary_data['blastp_result_count'] += 1
    cluster_summary[cluster_name] = summary_data

# --- Prepare and populate output DataFrame ---
# This section calculates percentages and other aggregated metrics for each cluster
# based on the summarized data, then constructs the final output DataFrame.
output_data = []
for cluster_name, summary in cluster_summary.items():
    total_count = summary['total_count']
    reviewed_count = summary['reviewed_count']

    reviewed_name_r = find_highest_value_item(summary['reviewed_name'], plddt_dict)

    eggnog_description_percentage = summary['eggnog_annotation_count'] / total_count if total_count > 0 else 0
    eggnog_pfam_percentage = summary['eggnog_pfam_count'] / total_count if total_count > 0 else 0
    eggnog_pname_percentage = summary['eggnog_pname_count'] / total_count if total_count > 0 else 0
    blastp_percentage = summary['blastp_result_count'] / total_count if total_count > 0 else 0
    unreviewed_percentage = (1 - (reviewed_count / total_count)) * 100 if total_count > 0 else 0

    species_name_count_str = '; '.join(
        [f"{spe_name}: {', '.join(elements)}" for spe_name, elements in summary['spe_name'].items()]
    )
    spe_count = len(summary['spe_name'])

    output_data.append({
        'Cluster': cluster_name,
        'species_name_count': species_name_count_str,
        'spe_count': spe_count,
        'cluster_size': total_count,
        #'reviewed_protein': reviewed_name_r,
        'unreviewed_percentage': unreviewed_percentage,
        #'eggnog_description_percentage': eggnog_description_percentage,
        #'eggnog_pfam_percentage': eggnog_pfam_percentage,
        #'eggnog_pname_percentage': eggnog_pname_percentage,
        #'blastp_result_percentage': blastp_percentage
    })
output_df = pd.DataFrame(output_data)


# --- Extract structure names and calculate pLDDT metrics ---
# This segment extracts specific structure names (especially for Arabidopsis) from
# the 'species_name_count' column, then finds the representative structure
# with the highest pLDDT value and calculates average pLDDT for all structures in a cluster.
def extract_structure_names_pattern(species_name_count_str, target_species=None):
    names = []
    pattern = r'AF-(.*?)-F'
    if not isinstance(species_name_count_str, str):
        return [], [] if target_species else []

    for item in species_name_count_str.split(';'):
        parts = item.split(':')
        if len(parts) == 2:
            current_spe_name = parts[0].strip()
            elements_str = parts[1].strip()
            for filename in elements_str.split(', '):
                match = re.search(pattern, filename)
                if match:
                    if target_species == current_spe_name:
                        names.append(match.group(1))
                    elif target_species is None:
                        names.append(match.group(1))
    return names


represent_list = []
represent_plddt_list = []
plddt_average = []
represent_eggnog_pfam_list = []
represent_eggnog_description_list = []
represent_eggnog_pname_list = []

for _, row in output_df.iterrows():
    all_structure_names = extract_structure_names_pattern(row['species_name_count'])
    arab_names = extract_structure_names_pattern(row['species_name_count'], target_species='UP000006548')

    max_plddt = float('-inf')
    max_structure = None
    scores_for_average = [plddt_dict[name] for name in all_structure_names if name in plddt_dict]

    if arab_names:
        for name in arab_names:
            if name in plddt_dict and plddt_dict[name] > max_plddt:
                max_plddt = plddt_dict[name]
                max_structure = name
    elif all_structure_names:
        for name in all_structure_names:
            if name in plddt_dict and plddt_dict[name] > max_plddt:
                max_plddt = plddt_dict[name]
                max_structure = name

    represent_list.append(max_structure if max_structure != float('-inf') else 'N/A')
    represent_plddt_list.append(max_plddt if max_plddt != float('-inf') else 'N/A')
    plddt_average.append(sum(scores_for_average) / len(scores_for_average) if scores_for_average else 'N/A')

    current_represent_pfam = 'N/A'
    current_represent_description = 'N/A'
    current_represent_pname = 'N/A'

    if max_structure and max_structure != 'N/A' and max_structure in entry_info_egg:
        egg_data = entry_info_egg[max_structure]
        temp_description = str(egg_data.get('description', ''))
        temp_pfam = str(egg_data.get('pfam', ''))
        temp_pname = str(egg_data.get('preferred_name', ''))

        if not exclude_pattern.search(temp_description) and temp_description != '-':
            current_represent_description = temp_description
        if not exclude_pattern.search(temp_pfam) and temp_pfam != '-':
            current_represent_pfam = temp_pfam
        if not exclude_pattern.search(temp_pname) and temp_pname != '-':
            current_represent_pname = temp_pname

    represent_eggnog_pfam_list.append(current_represent_pfam)
    represent_eggnog_description_list.append(current_represent_description)
    represent_eggnog_pname_list.append(current_represent_pname)

output_df['represent'] = represent_list
output_df['represent_plddt'] = represent_plddt_list
output_df['plddt_average'] = plddt_average
output_df['represent_eggnog_pfam'] = represent_eggnog_pfam_list
output_df['represent_eggnog_description'] = represent_eggnog_description_list
output_df['represent_eggnog_pname'] = represent_eggnog_pname_list

# --- Map UniProt IDs to AGI IDs and identify matching TAIR annotation files ---
# This final segment links representative UniProt IDs to their corresponding AGI IDs
# and checks which TAIR annotation files (based on common annotation categories)
# contain those AGI IDs, then stores this information in the results DataFrame.
uni_to_agi = {}
with open(uni2tair_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            uni_id = parts[0]
            agi_ids_str = parts[1]
            if ';' in agi_ids_str:
                uni_to_agi[uni_id] = [i.split('.')[0] for i in agi_ids_str.split(';')]
            else:
                uni_to_agi[uni_id] = [agi_ids_str.split('.')[0]]

tair_annotation_files = [f for f in os.listdir(tair_annotate_dir) if f.endswith('.txt') and 'AGI' not in f]
tair_file_contents = {}
for filename in tair_annotation_files:
    file_path = os.path.join(tair_annotate_dir, filename)
    with open(file_path, 'r') as f:
        tair_file_contents[filename.split('.')[0]] = f.read()

result_df = output_df.copy()
result_df['matching_files'] = None

for index, row in output_df.iterrows():
    uni_id = row['represent']
    agi_ids = uni_to_agi.get(uni_id, [])
    matching_tair_files = []

    if agi_ids:
        for agi_id in agi_ids:
            for file_base_name, content in tair_file_contents.items():
                if agi_id in content:
                    matching_tair_files.append(file_base_name)
        if matching_tair_files:
            result_df.at[index, 'matching_files'] = ', '.join(sorted(list(set(matching_tair_files))))

print(result_df)
result_df.to_csv('Stastic/Nonsingleton_clusters_allinfo_stastic.csv', index=False)
