# BeatAML2 Data
Text from the website: [BeatAML2 Github](https://biodev.github.io/BeatAML2/)
> Acute myeloid leukemia (AML) is a cancer of myeloid-lineage cells with limited therapeutic options. We previously 
> combined ex vivo drug sensitivity with genomic, transcriptomic, and clinical annotations for a large cohort of AML 
> patients, which facilitated discovery of functional genomic correlates. Here, we present a new dataset, which has been
> harmonized with our initial report to yield a cumulative cohort of 805 patients (942 specimens). We show strong 
> cross-cohort validation and identify new features of drug response. Further, deconvoluting transcriptomic data shows 
> that drug sensitivity is governed broadly by AML cell differentiation state, sometimes conditionally impacting other 
> correlates of response. Finally, modeling of clinical outcome reveals a single gene, PEAR1, to be among the strongest 
> predictors of patient survival, especially for young patients. Collectively, this report expands a large functional 
> genomic resource, offers new avenues for mechanistic exploration and drug development, and reveals new tools for 
> predicting outcome in AML.

Citation: Bottomly, D., Long, N., Schultz, A. R., Kurtz, S. E., Tognon, C. E., Johnson, K., … & Tyner, J. W. (2022). Integrative analysis of drug response and clinical outcome in acute myeloid leukemia. Cancer Cell, 40(8), 850-864.

-> They provide drug response data, gene expression data (RNAseq), and mutation data (whole exome sequencing).

## Response

Raw unfitted drug response values:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt

Original IC50 values from the study:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_probit_curve_fits_v4_dbgap.txt

They were reprocessed with [beataml_preprocessing.ipynb](response/beataml_preprocessing.ipynb).

## Other modalities
The other files were processed with [process_beataml_omics.ipynb](process_beataml_omics.ipynb).
The mapping file between all 3 modalities is stored in https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_sample_mapping.xlsx


### Gene expression

Gene Expression Values:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt

### Mutation

Mutation Values: 
https://github.com/biodev/beataml2.0_data/raw/main/beataml_wes_wv1to4_mutations_dbgap.txt

The processing steps are in the notebook but they have just measured 2294 genes. These genes do not have enough
overlap with our gene lists. We therefore chose not to include this data:
```
Missing landmark genes: 760/867
Missing drug target genes: 196/242
Missing paccman genes: 1701/1957
```

### Drug features
See [drugs.ipynb](drugs.ipynb) for drug_smiles.csv and drug fingerprints.
For the drug graphs and ChemBERTa embeddings go to drevalpy/datasets/featurizer/: 
```{bash}
python create_drug_graphs.py --data_path path_to_this_project BeatAML2
```
For the ChemBERTa embeddings: 
```{bash}
python create_drug_embeddings.py --data_path path_to_this_project BeatAML2
```
For MolGNet: 


