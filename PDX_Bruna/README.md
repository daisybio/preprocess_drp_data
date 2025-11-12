# PDX data by Bruna et al.

The PDX data was retrieved from the supplementary materials of the original publication:

Bruna, A., Rueda, O. M., Greenwood, W., Batra, A. S., Callari, M., Batra, R. N., ... & Caldas, C. (2016). A biobank of breast cancer explants with preserved intra-tumor heterogeneity to screen anticancer compounds. Cell, 167(1), 260-274.
DOI: 10.1016/j.cell.2016.08.041 

The data consists of breast cancer patient-derived tumor xenografts (PDTXs). 
From the PDTXs, short term cultures of cancer cells were established (PDTCs) and used for high-throughput drug screening:
"A selection of 22 different PDTX models were plated as PDTCs and 24 hr later screened with 108 compounds, representing a total of 6,634 drug tests performed [...] The effect of drug treatment on cell viability was determined by CellTiter-Glo".

"PDTCs derived from ten of the PDTX models were extensively characterized using WES, sWGS, and RNAexp (microarrays)."

"The PDTX samples were comprehensively molecularly characterized at several passages using sWGS (for CNAs), WES (for single nucleotide variations [SNVs]), reduced-representation bisulfite sequencing (“RRBS”) (for DNA methylation), and RNAexp (for global expression and pathway activity profiling)."

All files can be retrieved from https://figshare.com/s/4a3f6bc543e5ba85834c. 

Except for the drug responses, we use the processed data provided in the figshare repository.
The omics files were postprocessed with [process_bruna_data.ipynb](process_bruna_data.ipynb) to generate the final files used for modeling currently in the Zenodo.

## Response
The drug response data is provided in the file `RawDataDrugsSingleAgents.txt`.
It was reprocessed with CurveCurator in the notebook [Bruna_preprocess_raw_and_run_curvecurator.ipynb](response/Bruna_preprocess_raw_and_run_curvecurator.ipynb).

## Expression
The microarray file downloaded is `ExpressionSamples.txt`. The raw files would be available at EGAS00001001913 
but we chose to only re-process RNAseq data, so we worked with the provided values.

## Copy number variations
The corresponding file is `CNASamples.txt`. We are reprocessing it with [GISTIC2.0](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-4-r41) using the [GenePattern server](https://cloud.genepattern.org/gp/pages/index.jsf). 
The preparing code for GISTIC is in the Jupyter notebook.

We used the Human Hg19 reference genome and uploaded the segmentation file and used the default settings. The results were 
written to all_thresholded.by_genes.txt.
The remaining preprocessing code for Zenodo is in the Jupyter notebook, too.

## Mutation
The corresponding file is `SNVsSamples.txt`. We converted the file to boolean values (code in notebook).

## Methylation
The corresponding file is `PromoterMethylationSamples.txt` (RRBS).
It is just promoter methylation, so we do not have the locations on the chromosomes but gene names. One could, theoretically, 
map this using a GTF file. But ultimately, we decided not to include the Methylation data in the Zenodo, because the values
are not comparable (CpG vs. promoter methylation). Might do this at a later point in time. 

## Cell line names, drug names
See notebook, just extracted. For gene lists, see gene lists directory

## Drug features
See [drugs.ipynb](drugs.ipynb) for drug_smiles.csv and drug fingerprints.
For the drug graphs and ChemBERTa embeddings go to drevalpy/datasets/featurizer/: 
```{bash}
python create_drug_graphs.py --data_path path_to_this_project PDX_Bruna
```
For the ChemBERTa embeddings: 
```{bash}
python create_drug_embeddings.py --data_path path_to_this_project PDX_Bruna
```
For MolGNet: 

