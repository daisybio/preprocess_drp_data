import pandas as pd
import pubchempy as pcp

path_to_data = "/Users/judithbernett/PycharmProjects/drp_model_suite/data/"

datasets = {
    #"CCLE_curvecurator": f"{path_to_data}/CCLE/CCLE.csv",
    #        "GDSC1_curvecurator": f"{path_to_data}/GDSC1/GDSC1.csv",
    #        "GDSC2_curvecurator": f"{path_to_data}/GDSC2/GDSC2.csv",
    #        "CTRPv1_curvecurator": f"{path_to_data}/CTRPv1/CTRPv1.csv",
            "CTRPv2_curvecurator": f"{path_to_data}/CTRPv2/CTRPv2.csv",
            "Toy_Data": f"{path_to_data}/Toy_Data/toy_data.csv",
            }
results_df = pd.DataFrame(columns=["pubchem_id", "drug_name", "canonical_smiles", "cactvs_fingerprint", "fingerprint", "dataset"])
unmapped_drugs = dict()
for dataset, path in datasets.items():
    print(f"Processing {dataset}")
    tmp = pd.read_csv(path)
    tmp = tmp[['drug_name', 'pubchem_id']]
    tmp = tmp.drop_duplicates()
    tmp = tmp.reset_index(drop=True)
    #if dataset == "Toy_Data":
    #    toy_data_ids = dict()
    for idx, row in tmp.iterrows():
        if idx % 10 == 0:
            print(f"Processing drug {idx} out of {len(tmp)}")
        try:
            cid = row["pubchem_id"]
            drug_name = row["drug_name"]
            compound = pcp.get_compounds(cid, 'cid')[0]
            smiles = compound.canonical_smiles
            cactvs = compound.cactvs_fingerprint
            fingerprint = compound.fingerprint
            results_df = pd.concat([results_df, pd.DataFrame.from_dict({
                "pubchem_id": [cid],
                "drug_name": [drug_name],
                "canonical_smiles": [smiles],
                "cactvs_fingerprint": [cactvs],
                "fingerprint": [fingerprint],
                "dataset": [dataset]
            })])
            #if dataset == "Toy_Data":
            #    toy_data_ids[name] = drug_id
        except:
            print(f"Could not find CID for {row['drug_name']}")
            if row['drug_name'] not in unmapped_drugs:
                unmapped_drugs[row['drug_name']] = [dataset]
            else:
                unmapped_drugs[row['drug_name']].append(dataset)
    #if dataset == "Toy_Data":
    #    tmp["drug_id"] = [toy_data_ids[name] for name in tmp["drug_name"]]
    #    tmp.to_csv(f"{path_to_data}/Toy_Data/toy_data.csv", index=False)

results_df.to_csv("pubchem_id_mapping.csv", index=False)
unmapped_drugs_df = pd.DataFrame.from_dict(unmapped_drugs, orient="index")
unmapped_drugs_df.to_csv("unmapped_drugs.csv", index=True)