import pandas as pd
import pubchempy as pcp
import os

path_to_data = "/Users/judithbernett/PycharmProjects/drp_model_suite/data/"
"""
filename = "drug_name_to_demorgan_64_map.csv"

fingerprint_cid_mapping = dict()
fingerprints = pd.read_csv(f"{path_to_data}/CCLE/drug_fingerprints/{filename}")
fingerprints = fingerprints.drop(columns=["Unnamed: 0"])
fingerprint_drugs = set(fingerprints.columns)

for drug in fingerprint_drugs:
    try:
        compound = pcp.get_compounds(drug, 'name')[0]
        drug_id = compound.cid
        fingerprint_cid_mapping[drug] = drug_id
    except:
        print(f"Could not find CID for {drug}")

fingerprint_cid_mapping["Picolinici-acid"] = 1018
fingerprint_cid_mapping["Nutlin-3a (-)"] = 11433190
fingerprint_cid_mapping["KRAS (G12C) Inhibitor-12"] = 135148916
fingerprint_cid_mapping["MPS-1-IN-1"] = 25195352
fingerprint_cid_mapping["KIN001-260"] = 10451420
fingerprint_cid_mapping["Cetuximab"] = 85668777
fingerprint_cid_mapping["Venotoclax"] = 49846579
fingerprint_cid_mapping["JNK-9L"] = 25222038

fingerprint_cid_mapping_df = pd.DataFrame.from_dict(fingerprint_cid_mapping, orient="index", columns=["CID"])
fingerprint_cid_mapping_df.to_csv("drug_name_to_cid_map.csv")

# some are duplicated
duplicate_ids = fingerprint_cid_mapping_df[fingerprint_cid_mapping_df["CID"].isin(fingerprint_cid_mapping_df[fingerprint_cid_mapping_df.duplicated(subset=["CID"])]["CID"].tolist())]
print(duplicate_ids)
remove_duplicates = ["AZD8931", "AZD5363", "Venotoclax", "AZD6094"]
fingerprints = fingerprints.drop(columns=remove_duplicates)
# rename columns to CIDs
fingerprints = fingerprints.rename(columns=fingerprint_cid_mapping)
fingerprints.to_csv(f"{path_to_data}/CCLE/drug_fingerprints/{filename}", index=False)
"""
fingerprint_df = pd.read_csv("drug_name_to_cid_map.csv")
fingerprint_dict = fingerprint_df.set_index("Unnamed: 0").to_dict()
for filename in ["drug_name_to_demorgan_128_map.csv", "drug_name_to_demorgan_256_map.csv",
                 "drug_name_to_demorgan_512_map.csv", "drug_name_to_demorgan_1024_map.csv", "drug_name_to_demorgan_2048_map.csv"]:
    print(f"Mapping {filename}")
    fingerprints = pd.read_csv(f"{path_to_data}/CCLE/drug_fingerprints/{filename}")
    fingerprints = fingerprints.drop(columns=["Unnamed: 0"])
    remove_duplicates = ["AZD8931", "AZD5363", "Venotoclax", "AZD6094"]
    fingerprints = fingerprints.drop(columns=remove_duplicates)
    fingerprints = fingerprints.rename(columns=fingerprint_dict["CID"])
    fingerprints.to_csv(f"{path_to_data}/CCLE/drug_fingerprints/{filename}", index=False)
"""
dirname = "CCLE/DIPK_features/Drugs"
# list all filenames in the dir
filenames = os.listdir(f"{path_to_data}{dirname}")
fingerprint_cid_mapping = dict()
for file in filenames:
    drug_name = file.split(".")[0].split("_")[1]
    try:
        compound = pcp.get_compounds(drug_name, 'name')[0]
        drug_id = compound.cid
    except:
        manual_mappings = {"KIN001-260": 10451420,
                           "JNK-9L": 25222038,
                           "Nutlin-3a (-)": 11433190,
                           "Lestauritinib": 126565,
                           "MPS-1-IN-1": 25195352,
                           "Cetuximab": 85668777
                           }
        drug_id = manual_mappings[drug_name]
    fingerprint_cid_mapping[drug_name] = drug_id
    new_filename = f"{path_to_data}{dirname}/MolGNet_{drug_id}.csv"
    os.rename(f"{path_to_data}{dirname}/{file}", new_filename)


fingerprint_cid_mapping_df = pd.DataFrame.from_dict(fingerprint_cid_mapping, orient="index", columns=["CID"])
fingerprint_cid_mapping_df.to_csv("drug_name_to_cid_map_molgnet.csv")
"""