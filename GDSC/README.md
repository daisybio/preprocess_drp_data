# GDSC

1. [Transcriptomics](#transcriptomics)
2. [Mutations](#mutations)
3. [Methylation](#methylation)
4. [Copy number variation](#copy-number-variation)
5. [Proteomics](#proteomics)
6. [Response](#response)

## Transcriptomics
From the [GDSC Data Portal](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html),
the gene expression data (RMA normalised expression data for cell-lines, Cell_line_RMA_proc_basalExp.txt) and two 
annotations (methSampleId_2_cosmicIds.xlsx, Mapping between cell-line COSMIC identifiers and cell-line methylation data identifiers + 
TableS1E.xlsx, Annotation of cell lines used in the GDSC dataset) was downloaded. 

The data was transformed with the following code:
```{python}
import pandas as pd
gdsc_gex = pd.read_csv('GDSC/gene_expression/Cell_line_RMA_proc_basalExp.txt', sep='\t')
cosmic_ids = gdsc_gex.columns[2:].to_list()
cosmic_ids = [int(cosmic_id.split('.')[1]) for cosmic_id in cosmic_ids]

annotation = pd.read_excel('GDSC/annotation/methSampleId_2_cosmicIds.xlsx')
# drop columns where 'cosmic_id' is NaN
annotation = annotation.dropna(subset=['cosmic_id'])
anno_dict = {int(row['cosmic_id']): row['Sample_Name'] for index, row in annotation.iterrows()}

annotation2 = pd.read_excel('GDSC/annotation/TableS1E.xlsx', skiprows=[0,1,2,1005])
# drop first column which is empty
annotation2 = annotation2.drop(columns=annotation2.columns[0])
# read only the third row of the excel file
colnames = pd.read_excel('GDSC/annotation/TableS1E.xlsx', nrows=2, skiprows=1)
colnames = colnames.drop(columns=colnames.columns[0])
colnames = colnames.ffill()
colnames = colnames.iloc[1].tolist()
# replace '\n' in the column names with ' '
colnames = [col.replace('\n', ' ') for col in colnames]
annotation2.columns = colnames
update_dict = {int(row['COSMIC identifier']): row['Sample Name'] for index, row in annotation2.iterrows()}

anno_dict.update(update_dict)

unmapped_dict = {
    1298355: 'NCI-H2286',
    1723793: 'M000921',
    1723794: 'M980513',
    1240156: 'Ishikawa (Heraklio) 02 ER-',
    1298154: 'JHU-028',
    1659787: 'MET-2B'
}
anno_dict.update(unmapped_dict)

mapped_cosmic_ids = [str(anno_dict.get(cosmic_id, f'Cosmic:{cosmic_id}')) for cosmic_id in cosmic_ids]
gdsc_gex.columns = ['GENE', 'GENE_SYMBOL'] + mapped_cosmic_ids
gdsc_gex = gdsc_gex.drop(columns=['GENE_SYMBOL'])
gdsc_gex = gdsc_gex.set_index('GENE')
gdsc_gex = gdsc_gex.T
gdsc_gex.index.name = 'CELL_LINE_NAME'
gdsc_gex.to_csv('GDSC/gene_expression/reprocessed_gdsc_gex.csv')
```

Then, the data was mapped to cellosaurus IDs with the code in utils/convert_to_cello.py

## Mutations

The data was downloaded from [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads):
Mutation Data -> Mutations All. 

The data was filtered for **only coding** mutations which are **not silent** and then transformed with the following code:
```{python} 
import pandas as pd
mut = pd.read_csv('SangerCellModelPassports/mutation/mutations_all_20230202.csv')
# filter for coding mutations 
mut = mut[mut['coding'] == True]
# filter out silent mutations
mut = mut[mut['effect'] != 'silent']
# filter for 'gene_symbol' and 'model_id', introduce a new column: 'mutated' which is True everywhere
mut = mut[['gene_symbol', 'model_id']]
mut['mutated'] = True
# mapping for sanger
cellosaurus = pd.read_csv("mapping/cellosaurus_01_2024.csv")
cellosaurus = cellosaurus.fillna("")
cellosaurus_sid_dict = {}
cellosaurus_id_dict = {}
for index, row in cellosaurus.iterrows():
    organism = row["OX"]
    if 'Human' not in organism:
        continue
    if "Cell_Model_Passport" in row["DR"]:
        # split DR column by ',', iterate until you encounter an ID starting with SIDM
        for element in row["DR"].split(","):
            if element.startswith("SIDM"):
                cellosaurus_sid_dict[element] = row["AC"]
                cellosaurus_id_dict[element] = row["ID"]
    else:
        continue

mut['CELL_LINE_NAME'] = [cellosaurus_id_dict.get(model_id) for model_id in mut['model_id']]
mut['cellosaurus_id'] = [cellosaurus_sid_dict.get(model_id) for model_id in mut['model_id']]
mut = mut.drop(columns=['model_id'])
mut = mut.dropna(subset=['CELL_LINE_NAME'])
mut = mut.drop_duplicates()
mut = mut.set_index(["cellosaurus_id", "CELL_LINE_NAME"])
# pivot the table
mut = mut.pivot(columns='gene_symbol', values='mutated')
mut = mut.fillna(False)
mut = mut.reset_index()
mut.to_csv('SangerCellModelPassports/mutation/mutations_cellosaurus.csv', index=False)
```

## Methylation

From the [GDSC Data Portal](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html),
the methylation data (F2_METH_CELL_DATA, Pre-Processed beta values for all CpG islands across all the cell-line).

The data was transformed with the following code: 
```{python}
import pandas as pd
met = pd.read_csv('GDSC/methylation/F2_METH_CELL_DATA.txt', sep='\t')
met = met.set_index("Unnamed: 0")
met = met.T

annotation = pd.read_excel('GDSC/annotation/methSampleId_2_cosmicIds.xlsx')
anno_dict = {f'{row["Sentrix_ID"]}_{row["Sentrix_Position"]}': row["Sample_Name"] for index, row in annotation.iterrows()}
met["CELL_LINE_NAME"] = [anno_dict.get(id) for id in met.index]
met = met.set_index("CELL_LINE_NAME")
met.to_csv('GDSC/methylation/reprocessed_gdsc_met.csv')
```
Then, the data was mapped to cellosaurus IDs with the code in utils/convert_to_cello.py

## Copy number variation

The data was downloaded from [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads): 
Copy Number Data -> Copy Number (SNP6) -> PICNIC absolute copy numbers and GISTIC scores derived from Affymetrix SNP6.0 array data.

The data was transformed with the following code:
```{python}
import pandas as pd
cnv_sanger = pd.read_csv('SangerCellModelPassports/cnv/cnv_gistic_20191101.csv', skiprows=[0,2])
cnv_sanger = cnv_sanger.drop(columns=['model_name'])
cnv_sanger = cnv_sanger.set_index("Unnamed: 1")
cnv_sanger = cnv_sanger.T
cnv_sanger.index.name = 'CELL_LINE_NAME'
cnv_sanger = cnv_sanger.reset_index()
cnv_sanger.to_csv('GDSC/cnv/copy_number_variation_gistic.csv')
```

Then, it was mapped to cellosaurus IDs with the code in utils/convert_to_cello.py

## Proteomics

A DIA screen was done by [Gonçalves et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1535610822002744).
Their raw data is located in the [PRIDE: PXD030304](https://www.ebi.ac.uk/pride/archive/projects/PXD030304).
It contains 949 cell lines with 3-10 replicates. 

From this link, the table `ProCan-DepMapSanger_DIANN_output.tsv`(221 Gb) was downloaded. 

It was first filtered for FDR: 
```{python}
import pandas as pd
chunks = pd.read_csv('ProCan-DepMapSanger_DIANN_output.tsv', chunksize=1000000,sep='\t')
for i,df in enumerate(chunks):
    df2 = df[(df['Global.Q.Value'] <= 0.01) & (df['Proteotypic'] == 1)]
    if len(df2):
        df2.to_csv(f'./chunks/df{i}.tsv',index=False,sep='\t')
```
These chunks were concatenated:
```{bash}
head -n 1  df0.tsv > ../q_1_percent.tsv
for file in $(ls *.tsv)
do
                echo $file
                cat $file|tail -n +2 >>../q_1_percent.tsv
done
```
Then, the max. LFQ was calculated with the diann R package
```{R}
install.packages("devtools")
library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)

df <- diann_load('q_1_percent.tsv')
protein.norm <- diann_maxlfq(df, group.header="Protein.Ids", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
write.csv(protein.norm, "protein_matrix_maxlfq_diann-normalised.tsv")
```
The header columns rawfiles were mapped to their corresponding cell lines annotations with
```{python}
import pandas as pd
normalized_file_df = pd.read_csv('protein_matrix_maxlfq_diann-normalised.csv')
normalized_file_df.columns = normalized_file_df.columns.str.replace('./','')
normalized_file_df.columns = normalized_file_df.columns.str.replace('.dia','')
mapping_df = pd.read_excel('mapping_Data.xlsx',sheet_name='Replicate level sample info',skiprows=1)
mapping_df_dic = mapping_df[['Automatic_MS_filename','Cell_line']].set_index('Automatic_MS_filename')['Cell_line'].to_dict()
normalized_file_df.columns = normalized_file_df.columns.map(mapping_df_dic)
normalized_file_df.to_csv('mapped_protein_matrix_maxlfq_diann-normalised.csv',index=False)
```
Finally, iBAQ was calculated with
```{python}
import pandas as pd
ibaq_mapping_df = pd.read_csv('ibaq_mapping_df.csv')
max_lfq_df = pd.read_csv('protein_matrix_maxlfq_diann-normalised.csv')
max_lfq_df = max_lfq_df.set_index('Unnamed: 0')
merged_df = max_lfq_df.merge(ibaq_mapping_df[['uniprot_id',     'max_num_peptide']],right_on='uniprot_id',left_index=True)
merged_df = merged_df.loc[:,~merged_df.columns.str.contains('control|Control|Unnamed:')]
columns_list = list(merged_df.columns[~merged_df.columns.str.contains('uniprot_id|max_num_peptide')])
max_lfq_df = merged_df[columns_list]
pseudo_ibaq_df = max_lfq_df.div(merged_df['max_num_peptide'], axis=0)
#
max_lfq_df['uniprot_id'] = merged_df['uniprot_id']
max_lfq_df['num_peptides'] = merged_df['max_num_peptide']
pseudo_ibaq_df['uniprot_id'] = merged_df['uniprot_id']
#
max_lfq_df.to_csv('max_lfq_df_for_PrDB.csv')
pseudo_ibaq_df.to_csv('ibaq_df_for_PrDB.csv')
```

and then mapped to cellosaurus IDs with the code in utils/convert_to_cello.py.

Alternatively, processed data can be downloaded from [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads): 
Proteomics Data (DIA-MS) (averaged intensities and averaged z-scores over the replicates).


## Response

There are two response datasets available from https://www.cancerrxgene.org/downloads/bulk_download:

* GDSC1: GDSC1_fitted_dose_response_27Oct23.xlsx
* GDSC2: GDSC2_fitted_dose_response_27Oct23.xlsx

Those were mapped to cellosaurus IDs with the code in utils/convert_to_cello.py. 

Additionally, the raw data was reprocessed with CurveCurator. The raw data can be obtained from the same website.

