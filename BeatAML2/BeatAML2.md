# BeatAML2 Data
Citation: Bottomly, D., Long, N., Schultz, A. R., Kurtz, S. E., Tognon, C. E., Johnson, K., … & Tyner, J. W. (2022). Integrative analysis of drug response and clinical outcome in acute myeloid leukemia. Cancer Cell, 40(8), 850-864.
## Response

The original response data was retrieved from the [BeatAML2 Github](https://biodev.github.io/BeatAML2/)  using:

raw unfitted drug response values:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt

Gene Expression Values:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt

Original IC50 values from the study:
https://github.com/biodev/beataml2.0_data/raw/main/beataml_probit_curve_fits_v4_dbgap.txt


# harmonize drug ids
We process the raw unfitted drug response values in beataml_preprocessing.ipynb.
We first try to find as many pubchem_ids as possible by mapping the drug names to other screens and via synonym search etc.
When a PubChem ID cannot be found, we use the original drug name.


# create raw dose response file for curve curator
We then merge these mappings back into the BeatAML2 response tables to harmonize all drug identifiers across datasets.
"BeatAML2_raw.csv". This can be used for DrEval with curve curator option.

# feature dataset
we then create all the needed feature datasets (gene expression, fingerprints, drug an cell line names)


