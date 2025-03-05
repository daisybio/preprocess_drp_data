import pandas as pd
import numpy as np
import os

path_to_data = "/Users/judithbernett/PycharmProjects/drp_model_suite/data/"
dataset = "CTRPv2"
cross_dataset = "GDSC2"

# set seed
np.random.seed(42)

# Create a toy dataset
response = pd.read_csv(path_to_data + f"{dataset}/{dataset}.csv", dtype={"pubchem_id": str})
response_cross = pd.read_csv(path_to_data + f"{cross_dataset}/{cross_dataset}.csv", dtype={"pubchem_id": str})

def sample_random_ids(dataset1, dataset2, n):
    n_common = int(0.8 * n)
    n_exclusive = int(0.1 * n)
    ids_1 = set(dataset1)
    ids_2 = set(dataset2)
    random_ids = np.random.choice(list(ids_1.intersection(ids_2)), size=n_common, replace=False)
    if len(ids_1 - ids_2) > 0:
        random_ids = np.concatenate([random_ids, np.random.choice(list(ids_1 - ids_2), size=n_exclusive, replace=False)])
    if len(ids_2 - ids_1) > 0:
        random_ids = np.concatenate([random_ids, np.random.choice(list(ids_2 - ids_1), size=n_exclusive, replace=False)])
    return random_ids


random_drug_ids = sample_random_ids(response["pubchem_id"], response_cross["pubchem_id"], 40)
random_cellosaurus_ids = sample_random_ids(response["cellosaurus_id"], response_cross["cellosaurus_id"], 100)

print("Writing response data ...")
# Reponse data
for name, df in {"TOYv1": response, "TOYv2": response_cross}.items():
    df = df[df["pubchem_id"].isin(random_drug_ids)]
    df = df[df["cellosaurus_id"].isin(random_cellosaurus_ids)]
    df.to_csv(path_to_data + f"{name}/{name}.csv", index=False)
    cell_line_names = df[["cellosaurus_id", "cell_line_name"]]
    cell_line_names = cell_line_names.drop_duplicates()
    cell_line_names.to_csv(path_to_data + f"{name}/cell_line_names.csv", index=False)
    drug_names = df[["pubchem_id", "drug_name"]]
    drug_names = drug_names.drop_duplicates()
    drug_names.to_csv(path_to_data + f"{name}/drug_names.csv", index=False)

print("Writing drug features ...")
# Drug features
for name, ds in {"TOYv1": dataset, "TOYv2": cross_dataset}.items():
    fingerprints = pd.read_csv(path_to_data + f"{ds}/drug_fingerprints/pubchem_id_to_demorgan_128_map.csv")
    fingerprint_columns = [str(col) for col in fingerprints.columns if col in random_drug_ids]
    fingerprints = fingerprints[fingerprint_columns]
    # set last column to NA
    fingerprints.iloc[:, -1] = np.nan
    os.makedirs(path_to_data + f"{name}/drug_fingerprints", exist_ok=True)
    fingerprints.to_csv(path_to_data + f"{name}/drug_fingerprints/pubchem_id_to_demorgan_128_map.csv", index=False)

    dipk_features = os.listdir(path_to_data + f"{ds}/DIPK_features/Drugs")
    dipk_features = [f for f in dipk_features if f.split("_")[1].split(".csv")[0] in random_drug_ids]
    os.makedirs(path_to_data + f"{name}/DIPK_features/Drugs", exist_ok=True)
    for f in dipk_features:
        df = pd.read_csv(path_to_data + f"{ds}/DIPK_features/Drugs/{f}")
        df.to_csv(path_to_data + f"{name}/DIPK_features/Drugs/{f}", index=False)

# Cell line features
# toy genes: union of gene lists
toy_genes = set()
for gene_list in ["drug_target_genes_all_drugs.csv", "gene_list_paccmann_network_prop.csv", "landmark_genes.csv"]:
    genes = pd.read_csv(path_to_data + f"{dataset}/gene_lists/{gene_list}")
    toy_genes.update(genes["Symbol"])
toy_genes = list(toy_genes)

toy_genes_proteomics = set()
for gene_list in ["drug_target_genes_all_drugs_proteomics.csv", "gene_list_paccmann_network_prop_proteomics.csv", "landmark_genes_proteomics.csv"]:
    genes = pd.read_csv(path_to_data + f"{dataset}/gene_lists/{gene_list}")
    toy_genes_proteomics.update(genes["Symbol"])
toy_genes_proteomics = list(toy_genes_proteomics)

print("Writing gene expression data ...")
gex = pd.read_csv(path_to_data + f"{dataset}/gene_expression.csv")
gex_cross = pd.read_csv(path_to_data + f"{cross_dataset}/gene_expression.csv")

for name, ds in {"TOYv1": gex, "TOYv2": gex_cross}.items():
    ds = ds[ds["cellosaurus_id"].isin(random_cellosaurus_ids)]
    ds = ds[["cellosaurus_id", "cell_line_name"] + list(toy_genes)]
    ds.to_csv(path_to_data + f"{name}/gene_expression.csv", index=False)

print("Writing methylation data ...")
methylation = pd.read_csv(path_to_data + f"{dataset}/methylation.csv")
methylation_cross = pd.read_csv(path_to_data + f"{cross_dataset}/methylation.csv")
toy_methylation_features = sample_random_ids(methylation.columns[2:], methylation_cross.columns[2:], 100)

for name, ds in {"TOYv1": methylation, "TOYv2": methylation_cross}.items():
    ds = ds[ds["cellosaurus_id"].isin(random_cellosaurus_ids)]
    cols = list(set(ds.columns[2:]) & set(toy_methylation_features))
    ds = ds[["cellosaurus_id", "cell_line_name"] + cols]
    ds.to_csv(path_to_data + f"{name}/methylation.csv", index=False)

print("Writing other omics ...")
# these features are the same
for datatype in ["copy_number_variation_gistic", "mutations", "proteomics"]:
    df = pd.read_csv(path_to_data + f"{dataset}/{datatype}.csv")
    df = df[df["cellosaurus_id"].isin(random_cellosaurus_ids)]
    if datatype == "proteomics":
        df = df[["cellosaurus_id", "cell_line_name"] + toy_genes_proteomics]
    else:
        df = df[["cellosaurus_id", "cell_line_name"] + toy_genes]
    df.to_csv(path_to_data + f"TOYv1/{datatype}.csv", index=False)
    df.to_csv(path_to_data + f"TOYv2/{datatype}.csv", index=False)

print("Writing gene intersection lists ...")
# gene lists
for datatype in ["copy_number_variation_gistic", "gene_expression", "methylation", "mutations", "proteomics"]:
    df1 = pd.read_csv(path_to_data + f"TOYv1/{datatype}.csv")
    df2 = pd.read_csv(path_to_data + f"TOYv2/{datatype}.csv")
    genes = list(set(df1.columns[2:]) & set(df2.columns[2:]))
    genes = pd.DataFrame({"Symbol": genes})
    os.makedirs(path_to_data + f"TOYv1/gene_lists", exist_ok=True)
    os.makedirs(path_to_data + f"TOYv2/gene_lists", exist_ok=True)
    genes.to_csv(path_to_data + f"TOYv1/gene_lists/{datatype}_intersection.csv", index=False)
    genes.to_csv(path_to_data + f"TOYv2/gene_lists/{datatype}_intersection.csv", index=False)
