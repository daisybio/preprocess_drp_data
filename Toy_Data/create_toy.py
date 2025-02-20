import pandas as pd

path_to_data = "/Users/judithbernett/PycharmProjects/drp_model_suite/data/"
pubchem_ids_toy = pd.read_csv("pubchem_ids_toy.csv")
pubchem_ids_toy["pubchem_id"] = pubchem_ids_toy["pubchem_id"].astype(str)
cvcl_ids_toy = pd.read_csv("cvcl_ids_toy.csv")
toy_genes = pd.read_csv("toy_genes.csv")
toy_methylation_features = pd.read_csv("toy_methylation_features.csv")

# Create a toy dataset
gdsc1_response = pd.read_csv(path_to_data + "GDSC1/GDSC1.csv", dtype={"pubchem_id": str})
gdsc1_response = gdsc1_response[gdsc1_response["pubchem_id"].isin(pubchem_ids_toy["pubchem_id"])]
gdsc1_response = gdsc1_response[gdsc1_response["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
gdsc1_response.to_csv(path_to_data + "Toy_Data/toy_data.csv", index=False)

cell_line_names = gdsc1_response[["cellosaurus_id", "cell_line_name"]]
cell_line_names = cell_line_names.drop_duplicates()
cell_line_names.to_csv(path_to_data + "Toy_Data/cell_line_names.csv", index=False)

drug_names = gdsc1_response[["pubchem_id", "drug_name"]]
drug_names = drug_names.drop_duplicates()
drug_names.to_csv(path_to_data + "Toy_Data/drug_names.csv", index=False)

fingerprints = pd.read_csv(path_to_data + "GDSC1/drug_fingerprints/pubchem_id_to_demorgan_128_map.csv")
fingerprint_columns = pubchem_ids_toy["pubchem_id"].tolist()
fingerprints = fingerprints[fingerprint_columns]
fingerprints.to_csv(path_to_data + "Toy_Data/drug_fingerprints/pubchem_id_to_demorgan_128_map.csv", index=False)

gdsc1_cnv = pd.read_csv(path_to_data + "GDSC1/copy_number_variation_gistic.csv")
gdsc1_cnv = gdsc1_cnv[gdsc1_cnv["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
all_columns = ["cellosaurus_id", "cell_line_name"] + toy_genes["gene_id"].tolist()
gdsc1_cnv = gdsc1_cnv[all_columns]
gdsc1_cnv.to_csv(path_to_data + "Toy_Data/copy_number_variation_gistic.csv", index=False)

gdsc1_gex = pd.read_csv(path_to_data + "GDSC1/gene_expression.csv")
gdsc1_gex = gdsc1_gex[gdsc1_gex["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
gdsc1_gex = gdsc1_gex[all_columns]
gdsc1_gex.to_csv(path_to_data + "Toy_Data/gene_expression.csv", index=False)

gdsc1_methylation = pd.read_csv(path_to_data + "GDSC1/methylation.csv")
gdsc1_methylation = gdsc1_methylation[gdsc1_methylation["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
all_columns_met = ["cellosaurus_id", "cell_line_name"] + toy_methylation_features["methylation_feature"].tolist()
gdsc1_methylation = gdsc1_methylation[all_columns_met]
gdsc1_methylation.to_csv(path_to_data + "Toy_Data/methylation.csv", index=False)

gdsc1_mutation = pd.read_csv(path_to_data + "GDSC1/mutations.csv")
gdsc1_mutation = gdsc1_mutation[gdsc1_mutation["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
mutation_features = set(gdsc1_mutation.columns) - {"cellosaurus_id", "cell_line_name"}
mutation_features = mutation_features.intersection(set(toy_genes["gene_id"].tolist()))
all_columns_mut = ["cellosaurus_id", "cell_line_name"] + list(mutation_features)
gdsc1_mutation = gdsc1_mutation[all_columns_mut]
gdsc1_mutation.to_csv(path_to_data + "Toy_Data/mutations.csv", index=False)

gdsc1_proteomics = pd.read_csv(path_to_data + "GDSC1/proteomics.csv")
gdsc1_proteomics = gdsc1_proteomics[gdsc1_proteomics["cellosaurus_id"].isin(cvcl_ids_toy["cellosaurus_id"])]
prot_features = set(gdsc1_proteomics.columns) - {"cellosaurus_id", "cell_line_name"}
prot_features = prot_features.intersection(set(toy_genes["gene_id"].tolist()))
all_columns_prot = ["cellosaurus_id", "cell_line_name"] + list(prot_features)
gdsc1_proteomics = gdsc1_proteomics[all_columns_prot]
gdsc1_proteomics.to_csv(path_to_data + "Toy_Data/proteomics.csv", index=False)
