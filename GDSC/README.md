# GDSC

## Gene expression
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


## Methylation

From the [GDSC Data Portal](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html),
the methylation data (F2_METH_CELL_DATA, Pre-Processed beta values for all CpG islands across all the cell-line).

## Copy number variation

The data was downloaded from [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads): 
Copy Number Data -> Copy Number (SNP6) -> PICNIC absolute copy numbers and GISTIC scores derived from Affymetrix SNP6.0 array data.

The data was transformed with the following code:
```{python}
import pandas as pd
cnv_sanger = pd.read_csv('SangerCellModelPassports/cnv/cnv_20191101/cnv_gistic_20191101.csv', skiprows=[0,2])
cnv_sanger = cnv_sanger.drop(columns=['model_name'])
cnv_sanger = cnv_sanger.set_index("Unnamed: 1")
cnv_sanger = cnv_sanger.T
cnv_sanger.index.name = 'CELL_LINE_NAME'
cnv_sanger = cnv_sanger.reset_index()
```

Then, it was mapped to cellosaurus IDs with the code in utils/convert_to_cello.py