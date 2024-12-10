import pandas as pd
import numpy as np

path_to_GDSC_methylation = '../../GDSC/methylation/methylation.csv'

# Define match_max function
def match_max(from_start, to_end, chromosome, comp_data):
    overlaps = comp_data[(comp_data['Chromosome'] == chromosome) &
                         (to_end >= comp_data['Start']) &
                         (from_start <= comp_data['End'])]
    if overlaps.empty:
        return None, None
    overlaps['overlap'] = np.minimum(to_end, overlaps['End']) - np.maximum(from_start, overlaps['Start'])
    max_overlap = overlaps.loc[overlaps['overlap'].idxmax()]
    return max_overlap['V1'], max_overlap['overlap']


# Load methylation data
methylation_cpg = pd.read_csv("CCLE_RRBS_TSS_CpG_clusters_20180614.txt",
                              sep="\t", na_values='    NaN')


# Define function to extract 'from-to' information
def extract_from_to(row):
    from_to = row.split(";")
    from_part = from_to[0]
    to_part = from_to[-1]
    chrom, from_start = from_part.split(":")
    _, to_end = to_part.split(":")
    return f"chr{chrom}:{from_start}-{to_end}"


# Extract 'from-to' columns
methylation_cpg['from_to'] = methylation_cpg['CpG_sites_hg19'].apply(extract_from_to)
methylation_cpg[['Chromosome', 'Rest']] = methylation_cpg['from_to'].str.split(":", expand=True)
methylation_cpg[['Start', 'End']] = methylation_cpg['Rest'].str.split("-", expand=True)
methylation_cpg['Start'] = pd.to_numeric(methylation_cpg['Start'])
methylation_cpg['End'] = pd.to_numeric(methylation_cpg['End'])
methylation_cpg['Length'] = methylation_cpg['End'] - methylation_cpg['Start']

# Filter based on Length and Chromosome
methylation_cpg = methylation_cpg[(methylation_cpg['Length'] > 200) & (~methylation_cpg['Chromosome'].isin(['chrX', 'chrY']))]
methylation_cpg = methylation_cpg.drop_duplicates()

# Load comparison methylation data
comp_methylation = pd.read_csv(path_to_GDSC_methylation)
from_to_comp = pd.DataFrame({'V1': comp_methylation.columns[2:]})
from_to_comp[['Chromosome', 'Rest']] = from_to_comp['V1'].str.split(":", expand=True)
from_to_comp[['Start', 'End']] = from_to_comp['Rest'].str.split("-", expand=True)
from_to_comp['Start'] = pd.to_numeric(from_to_comp['Start'])
from_to_comp['End'] = pd.to_numeric(from_to_comp['End'])
from_to_comp['Length'] = from_to_comp['End'] - from_to_comp['Start']

# Apply match_max function
results = methylation_cpg.apply(
    lambda row: match_max(row['Start'], row['End'], row['Chromosome'], from_to_comp),
    axis=1
)
methylation_cpg['best_match'], methylation_cpg['overlap'] = zip(*results)

# Filter and reorder columns
cols = methylation_cpg.columns.tolist()
cols = cols[-2:] + cols[:-2]
methylation_cpg = methylation_cpg[cols]

# drop rows with no match
new_met_df = methylation_cpg.dropna(subset=['best_match'])
# fill NaN values with 0.0
new_met_df.fillna(0.0, inplace=True)
new_met_df = new_met_df.drop(columns=['Rest', 'Chromosome', 'Start', 'End', 'cluster_id', 'RefSeq_id', 'CpG_sites_hg19'])
new_met_df['best_match'] = new_met_df['best_match'].astype(str)
new_met_df['overlap'] = new_met_df['overlap'].astype(float)

# Remove duplicates based on best_match
unique_ids = new_met_df['best_match'].unique()
final_rows = dict()

for dup_id in unique_ids:
    matching_row = new_met_df[new_met_df['best_match'] == dup_id]
    if len(matching_row) > 1:
        # get the row with the highest overlap
        matching_row = matching_row[matching_row['overlap'] == matching_row['overlap'].max()]
        # if there are still duplicates, get the one with the highest avg_coverage
        if len(matching_row) > 1:
            matching_row = matching_row[matching_row['avg_coverage'] == matching_row['avg_coverage'].max()]
        # if there are still duplicates, get the first one
    # convert to Series
    matching_row = matching_row.iloc[0]
    final_rows[dup_id] = matching_row.tolist()

new_met_df = pd.DataFrame.from_dict(final_rows, orient='index', columns=new_met_df.columns)

new_met_df = new_met_df.drop(columns=['overlap', 'gene_name', 'avg_coverage', 'best_match', 'from_to', 'Length'])
# sort by index
new_met_df = new_met_df.T
# reset index
new_met_df = new_met_df.reset_index(names='CELL_LINE')
new_met_df.to_csv('methylation_pre_mapping.csv', index=False)
