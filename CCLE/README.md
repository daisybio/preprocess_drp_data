# CCLE Data

1. [Transcriptomics](#transcriptomics)
2. [Mutations](#mutations)
3. [Methylation](#methylation)
4. [Copy number variation](#copy-number-variation-gistic20)
5. [Proteomics](#proteomics)
6. [Response](#response)
7. [Gene lists](#gene-lists)
8. [Drug fingerprints](#drug-fingerprints)

## Transcriptomics

### Processed data
The already processed data can be downloaded from [DepMap](https://depmap.org/portal/data_page/?tab=allData): 
TPMs, processed with RSEM. There is also another file with TPMs stored at DepMap.

Another possible source is [CellModelPassports](https://cellmodelpassports.sanger.ac.uk/downloads). 
There, the CCLE Data (Broad Institute) is stored as well as the data from the Sanger institute. 
This data has been processed with the iRAP pipeline.

Furthermore, there is raw Microarray data available at [GEO](https://www.google.com/url?q=https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc%3DGSE36139&sa=D&source=editors&ust=1733326672230985&usg=AOvVaw34q-kblgqEp_0fIOJ_xzRx).

### Re-processing the RNAseq data with nf-core/rnaseq to get gene-level TPMs
Because none of these matrices agree with each other, we have decided to **re-process the raw RNAseq data from scratch.**

The CCLE cell lines were profiles using RNA-seq by [Ghandi et al. in 2019](10.1038/s41586-019-1186-3). 
The data is available via the NCBI SRA Run selector as BioProject PRJNA523380. 
From here, all RNA-seq metadata and the SRR accession numbers for the runs were downloaded.

The data was downloaded using SRA tools: 
```{bash}
cd /path/to/utils/destination
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz

export PATH=$PATH:$PWD/sratoolkit.3.0.1-ubuntu64/bin
cd /path/to/data/destination

prefetch --option-file SRR_Acc_List_CCLE.txt
``` 
Then, the FASTQ files were generated and gzipped using this script:

```{bash}
#!/bin/bash
while read run
do
	echo $run
	fasterq-dump $run
	echo gzipping
	gzip $run*.fastq
done < SRR_Acc_List_CCLE.txt
```
Then, the directories were deleted.

Afterward, the FASTQ files were processed using the [nf-core/RNA-seq](https://nf-co.re/rnaseq/3.14.0/) pipeline using this command: 
```{bash}
nextflow run nf-core/rnaseq --input rnaseq/samplesheet.csv --outdir /path/to/data/destination/nf_core/ --multiqc_title CCLE_star_salmon -c /path/to/data/destination/nextflow.config -profile singularity,slurm --fasta /path/to/star/genome.fa --gtf /path/to/ensembl107_GRCh38/Homo_sapiens.GRCh38.107.gtf --star_index /path/to/star/index -r 3.10.1
```

This was run in batches of 200 because otherwise, the work directory got too large. 
The salmon summary files got overwritten, so I ran the nf-core scripts to produce the summary files 
([`salmon_tximport.r`](gene_expression/salmon_tximport.r) because the results for each file are in each subfolder. 
Here, I also provided the correct Cellosaurus cell line names instead of the SRA Run Identifiers.

For prediction, we want to use the gene-level TPMs in ccle_salmon.gene_tpm_cvcl.tsv. 
It was postprocessed with
```{python}
import pandas as pd

cellosaurus = pd.read_csv("mapping/cellosaurus_01_2024.csv")
cellosaurus = cellosaurus.fillna("")
id_to_ac_dict = dict(zip(cellosaurus["AC"], cellosaurus["ID"]))

gex = pd.read_csv("CCLE/gene_expression/ccle_salmon.gene_tpm_cvcl.tsv", sep="\t")
# map ensg to gene symbol
mapping = pd.read_csv("CCLE/gene_expression/ensg_to_name.csv")
mapping_dict = dict(zip(mapping["initial_alias"], mapping["name"]))
gex["gene_name"] = [mapping_dict.get(gene_id) for gene_id in gex["gene_id"]]
gex = gex.dropna(subset=["gene_name"])
gex = gex.set_index("gene_name")
gex = gex.drop(columns=["gene_id"])
gex = gex.T
# map cell line names to cellosaurus IDs
gex["cell_line_name"] = [id_to_ac_dict.get(cell_line_name) for cell_line_name in gex.index]
gex = gex.reset_index()
gex = gex.rename(columns={"index": "cellosaurus_id"})
gex = gex.set_index(["cellosaurus_id", "cell_line_name"])
gex.to_csv("CCLE/gene_expression/salmon_gene_tpms.csv")
```

This **salmon_gene_tpms.csv** file is currently the file in the Zenodo **CCLE/gene_expression.csv**.

## Mutations

The mutation data was downloaded from [DepMap](https://depmap.org/portal/data_page/?tab=allData): Release 22Q2, file: CCLE_mutations.csv. 
It was transformed into the appropriate format with [process_mutation.R](mutation/process_mutation.R) to exclude 
**intronic, silent**, and mutations in **5'/3'UTR**.

Alternatively, the SangerCellModelPassports data can be used. For this, see the GDSC section. 
This data is **currently in Zenodo for GDSC1, GDSC2, and CCLE**.

## Methylation

From [DepMap](https://depmap.org/portal/data_page/?tab=allData): Methylation (RRBS), the data for the promoter CpG clusters was
downloaded (CCLE_RRBS_TSS_CpG_clusters_20180614.txt). 

This was preprocessed with [preprocess_methylation.R](methylation/preprocess_methylation.R). Here, the CpG cluster regions 
were matched to the ones from GDSC such that they share their variables. This was done by finding the regions that overlap 
most between the two datasets.

This was postprocessed in convert_to_cello.py to map the cell line names to Cellosaurus IDs.
It is now available in the Zenodo as **CCLE/methylation.csv**.

## Copy number variation: GISTIC2.0

The data can be downloaded from [Sanger Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads):
Copy Number Data -> Copy Number (SNP6) -> PICNIC absolute copy numbers and GISTIC scores derived from Affymetrix SNP6.0 array data
(see GDSC part).
This data is **currently in Zenodo for GDSC1, GDSC2, and CCLE**.

Alternatively, the copy number variation data can be downloaded from [DepMap](https://depmap.org/portal/data_page/?tab=allData): 
CCLE_copynumber_2013-12-03.seg.txt. 

We want to process them with [GISTIC2.0](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-4-r41) using 
the [GenePattern server](https://cloud.genepattern.org/gp/pages/index.jsf). 

We used the Human Hg19 reference genome and uploaded the segmentation file and used the default settings. The results were 
written to all_thresholded.by_genes.txt.

The data was preprocessed with the following script:
```{python}
import pandas as pd
cnv = pd.read_csv('CCLE/cnv/all_thresholded.by_genes.txt', sep='\t')
cnv = cnv.set_index("Gene Symbol")
cnv = cnv.drop(columns=['Locus ID', 'Cytoband'])
cnv = cnv.T
cnv.index.name = 'CELL_LINE_NAME'
cnv = cnv.reset_index()
cnv['CELL_LINE_NAME'] = cnv['CELL_LINE_NAME'].str.split('_').str[0]
cnv.to_csv('CCLE/cnv/reprocessed_cnv.csv', index=False)
```

The cell line IDs were mapped to Cellosaurus IDs with the code in [utils/convert_to_cello.py](utils/convert_to_cello.py) to produce 
``cnv_cellosaurus.csv``.

### Overlap between GISTIC dataset from Sanger & reprocessed dataset from DepMap

* In the Sanger Dataset, CNV data was measured for **978** cell lines and **20,669** genes.
* From the DepMap dataset, we get **1,040** cell lines (**1,036 unique**) and **23,109** genes. 
* From those, only **639** cell lines and **19,276** genes are overlapping.

Comparison: 
* Gene-wise correlation is alright: **0.771±0.067**
* Out of the 637 cell lines, on average, only 235.85±21.5 have the same value (~37%) per gene.

## Proteomics

We can also use the DIA-MS proteomics data from SangerCellModelPassports ([Gonçalves et al. (2021)](https://www.sciencedirect.com/science/article/pii/S1535610822002744)) 
for CCLE (see GDSC section).

This data is **currently in Zenodo for GDSC1, GDSC2, and CCLE**.

Alternatively, there was a CCLE proteomics TMT experiment done by [Nusinow et al.](https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419313856%3Fshowall%3Dtrue).
It contains 375 cell lines with no replicates. 
The processed dataset is available at [DepMap](https://depmap.org/portal/data_page/?tab=allData) or at [https://gygi.med.harvard.edu/publications/ccle](https://gygi.med.harvard.edu/publications/ccle). The raw data is at [MassIVE MSV000085836](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=02cd1b6a7c674f3ebdbed300b5d9aa57).

We obtained it >#TODO somehow< and mapped it to cellosaurus IDs with the code in utils/convert_to_cello.py.

## Response

The viability data was downloaded from the Supplementary Material (NIHMS361223-supplement-4.xlsx) 
of [Barretina et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3320027/#S2).

The unmodified version was obtained with
```{python}
import pandas as pd
import numpy as np
raw_df = pd.read_excel('response/NIHMS361223-supplement-4.xlsx', sheet_name=11, skiprows=2)
raw_df['LN_IC50'] = np.log(raw_df['IC50 (µM)'])
raw_df.rename(columns={'Primary Cell Line Name': 'CELL_LINE_NAME', 'Compound': 'DRUG_NAME'}, inplace=True)
raw_df.to_csv('response/response_CCLE.csv', index=False)
```
This was mapped to Cellosaurus IDs with the code in [utils/convert_to_cello.py](utils/convert_to_cello.py) to produce 
``response_CCLE_cellosaurus.csv``. It was uploaded to Zenodo as `response_CCLE.csv`.

It was reprocessed with CurveCurator: [CCLE_preprocess_raw_and_run_curvecurator.ipynb](response%2FCCLE_preprocess_raw_and_run_curvecurator.ipynb).
The data was filtered such that all pEC50s lie in the measured range. The resulting file `CCLE.csv` was uploaded to Zenodo.

The CTRP screens are larger:
- CTRPv1: 
- CTRPv2: 907 cell lines, 545 drugs. 



## Gene lists

## Drug fingerprints

