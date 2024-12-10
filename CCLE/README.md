# CCLE Data

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

## Mutations

The mutation data was downloaded from [DepMap](https://depmap.org/portal/data_page/?tab=allData): Release 22Q2, file: CCLE_mutations.csv. 
It was transformed into the appropriate format with [process_mutation.R](mutation/process_mutation.R).

## Methylation
From [DepMap](https://depmap.org/portal/data_page/?tab=allData): Methylation (RRBS), the data for the promoter CpG clusters was
downloaded (CCLE_RRBS_TSS_CpG_clusters_20180614.txt). 

This was preprocessed with [preprocess_methylation.R](methylation/preprocess_methylation.R). Here, the CpG cluster regions 
were matched to the ones from GDSC such that they share their variables. This was done by finding the regions that overlap 
most between the two datasets.

## Copy number variation: GISTIC2.0

The copy number variation data can be downloaded from [DepMap](https://depmap.org/portal/data_page/?tab=allData): 
CCLE_ABSOLUTE_combined_20181227.xlsx, CCLE ABSOLUTE copy number analysis results. 

We want to process them with [GISTIC2.0](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-4-r41)

## Proteomics
