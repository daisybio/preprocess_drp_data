# CTRP data

The CTRPv1 and CTRPv2 scans were mostly done on the CCLE cell lines. Hence, all OMICs are the same as the ones used for CCLE.
Only the response data differs.

## Response

The original response data was retrieved from the CTD^2 data portal using:

* https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv1.0_2013_pub_Cell_154_1151/ for CTRPv1
* https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/ for CTRPv2

It was reprocessed with CurveCurator: 
- [CTRPv1_preprocess_raw_and_run_curvecurator.ipynb](response%2FCTRPv1_preprocess_raw_and_run_curvecurator.ipynb)
- [CTRPv2_preprocess_raw_and_run_curvecurator.ipynb](response%2FCTRPv2_preprocess_raw_and_run_curvecurator.ipynb)

The reprocessed data was uploaded to Zenodo as `CTRPv1.csv` and `CTRPv2.csv`.

## cell line names

For the Naive models, we need the cell line names as features.
```{python}
import pandas as pd
response_ctrp = pd.read_csv('response/CTRPv1.csv') # pd.read_csv('response/CTRPv2.csv')
response_ctrp = response_ctrp[['cell_line_id', 'cell_line_name']]
response_ctrp = response_ctrp.rename(columns={'cell_line_id': 'cellosaurus_id'})
response_ctrp = response_ctrp.drop_duplicates()
response_ctrp.to_csv('cell_line_names_v1.csv', index=False) # to_csv('cell_line_names_v2.csv', index=False)
```
**242** unique cell lines for CTRPv1 and **886** unique cell lines for CTRPv2.

## drug names

```{python}
import pandas as pd
response_ctrp = pd.read_csv('response/CTRPv1.csv') # pd.read_csv('response/CTRPv2.csv')
response_ctrp = response_ctrp[['drug_id', 'cpd_name']]
response_ctrp = response_ctrp.rename(columns={'cpd_name': 'drug_name'})
response_ctrp = response_ctrp.drop_duplicates()
response_ctrp.to_csv('drug_names_v1.csv', index=False) # to_csv('drug_names_v2.csv', index=False)
```
**354** unique drugs for CTRPv1 and **573** unique drugs for CTRPv2.

