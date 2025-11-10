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