import pandas as pd
import numpy as np
from difflib import SequenceMatcher
import warnings

set_of_unmatched_cell_lines = set()
matched_cell_lines = dict()
renamed_cell_lines = dict()
no_match = set()


def create_cl_dict(df):
    # iterate over cellosaurus and make a dictionary mapping the 'ID' column to the 'AC' column and every cell line name
    # in the 'SY' column to the 'AC' column
    print("Creating cellosaurus dictionary ...")
    cello_ac_to_id_dict = {}
    cellosaurus_ac_dict = {}
    cellosaurus_sy_dict = {}
    species_dict = {}
    for index, row in df.iterrows():
        cello_ac_to_id_dict[row["AC"]] = row["ID"]
        # add species to species_dict
        species_dict[row["AC"]] = row["OX"]
        # check whether cellosaurus_dict[row['ID']] already exists, if not add it, if yes print the cell line name
        if row["ID"] not in cellosaurus_ac_dict and row["ID"] != "":
            cellosaurus_ac_dict[row["ID"]] = row["AC"]
        for cell_line_name in row["SY"].split("; "):
            # check whether cellosaurus_dict[cell_line_name] already exists, if not add it, if yes print the cell
            # line name
            if cell_line_name not in cellosaurus_sy_dict and cell_line_name != "":
                cellosaurus_sy_dict[cell_line_name] = row["AC"]
    return cellosaurus_ac_dict, cellosaurus_sy_dict, species_dict, cello_ac_to_id_dict


def create_cl_dict_cell_passp(df):
    # iterate over cellosaurus and make a dictionary mapping the 'ID' column to the 'SID' IDs from the DR column if
    # it contains Cell_Model_Passport
    print("Creating cellosaurus dictionary ...")
    cellosaurus_sid_dict = {}
    species_dict = {}
    for index, row in df.iterrows():
        # add species to species_dict
        species_dict[row["AC"]] = row["OX"]
        # check whether the DR column contains Cell_Model_Passport
        if "Cell_Model_Passport" in row["DR"]:
            # split DR column by ',', iterate until you encounter an ID starting with SIDM
            for element in row["DR"].split(","):
                if element.startswith("SIDM"):
                    cellosaurus_sid_dict[element] = row["AC"]
        else:
            continue
    return cellosaurus_sid_dict, species_dict


def map_to_cellosaurus(df, cl_dict_ac, cl_dict_sy, species_dict, cello_ac_to_id_dict, output_path):
    # iterate over the cell line names in the dataframe. Try to get the cell line name from the cellosaurus
    # dictionary. If it exists, put it in a column called 'cellosaurus_id'. If it doesn't exist, print the cell line
    # name
    print("Mapping cell line names to cellosaurus IDs ...")
    for index, row in df.iterrows():
        cell_line_name = row["cell_line_name"]
        try:
            df.loc[index, "cellosaurus_id"] = cl_dict_ac[cell_line_name]
            matched_cell_lines[cell_line_name] = cl_dict_ac[cell_line_name]
            species = species_dict.get(cl_dict_ac[cell_line_name])
            # if 'Human' is not part of species string, warn
            if "Human" not in species:
                warnings.warn(f"Cell line {cl_dict_ac[cell_line_name]} matched to {cell_line_name} is not human, but {species}.")
        except KeyError:
            try:
                cello_id = cl_dict_sy[cell_line_name]
                ac = cello_ac_to_id_dict[cello_id]
                df.loc[index, "cellosaurus_id"] = cl_dict_sy[cell_line_name]
                # rename cell_line_name to the one in the cellosaurus dictionary
                df["cell_line_name"] = df["cell_line_name"].replace({cell_line_name: ac})
                matched_cell_lines[cell_line_name] = cl_dict_sy[cell_line_name]
                species = species_dict.get(cl_dict_sy[cell_line_name])
                # if 'Human' is not part of species string, warn
                if "Human" not in species:
                    warnings.warn(f"Cell line {cl_dict_sy[cell_line_name]} matched to {cell_line_name} is not human, but {species}.")
                if cell_line_name not in set_of_unmatched_cell_lines:
                    set_of_unmatched_cell_lines.add(cell_line_name)
                    if (
                        SequenceMatcher(
                            a=cell_line_name.lower(),
                            b=list(cl_dict_ac.keys())[list(cl_dict_ac.values()).index(cl_dict_sy[cell_line_name])].lower(),
                        ).ratio()
                        < 0.7
                    ):
                        print(
                            f"no main match for {cell_line_name}, matched it to {cl_dict_sy[cell_line_name]} = {list(cl_dict_ac.keys())[list(cl_dict_ac.values()).index(cl_dict_sy[cell_line_name])]}, "
                            f"but the similarity was rather low:"
                            f"{SequenceMatcher(a=cell_line_name, b=list(cl_dict_ac.keys())[list(cl_dict_ac.values()).index(cl_dict_sy[cell_line_name])]).ratio()}"
                        )
            except KeyError:
                print(f"no match at all for {cell_line_name}")
                df.loc[index, "cellosaurus_id"] = pd.NA
                no_match.add(cell_line_name)

    # drop all rows where no cellosaurus ID could be found
    df = df.dropna(subset=["cellosaurus_id"])
    # save the gene expression dataframe with the cellosaurus IDs
    print("Saving dataframe with cellosaurus IDs ...")
    # make index to column 'cell_line_name' and make 'cellosaurus_id' the index
    df = df.set_index("cellosaurus_id")
    df.to_csv(output_path)


def map_to_cellosaurus_model_passp(df, cl_dict_sid, species_dict, output_path, ignore_columns=[]):
    # iterate over the cell line names in the dataframe. Try to get the cell line name from the cellosaurus
    # dictionary. If it exists, put it in a column called 'cellosaurus_id'. If it doesn't exist, print the cell line
    # name
    print("Mapping cell line names to cellosaurus IDs ...")
    # get all column names except the ones in ignore_columns
    columns = [col for col in df.columns if col not in ignore_columns]
    new_columns = list(df.columns.values)
    for cl_name in columns:
        try:
            new_columns[new_columns.index(cl_name)] = cl_dict_sid[cl_name]
            matched_cell_lines[cl_name] = cl_dict_sid[cl_name]
            species = species_dict.get(cl_dict_sid[cl_name])
            # if 'Human' is not part of species string, warn
            if "Human" not in species:
                warnings.warn(f"Cell line {cl_dict_sid[cl_name]} matched to {cl_name} is not human, but {species}.")
        except KeyError:
            print(f"no match at all for {cl_name}")
            no_match.add(cl_name)

    # rename columns
    df.columns = new_columns
    # save the gene expression dataframe with the cellosaurus IDs
    print("Saving dataframe with cellosaurus IDs ...")
    # make first column index
    df = df.set_index(new_columns[0])
    df.to_csv(output_path)


def preprocess_df(df, cl_col_name, renamed_dict):
    df = df.rename(columns={cl_col_name: "cell_line_name"})
    df["cell_line_name"] = df["cell_line_name"].replace(renamed_dict)
    renamed_cell_lines.update(renamed_dict)
    return df


def collapse_ln_ic50s(values):
    # for each value: take exponential of value. Take mean of all values. Take log of mean.
    return np.log(np.mean(np.exp(values)))


def preprocess_gex_gdsc():
    # read in gene expression dataframe
    print("Preprocessing gene expression dataframe ...")
    gex = pd.read_csv("../GDSC/gene_expression/reprocessed_gdsc_gex.csv", index_col=0)
    renamed = {
        "786-0": "786-O",
        "C32": "C32 [Human melanoma]",
        "G-292 Clone A141B1": "G-292 clone A141B1",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HH": "HH [Human lymphoma]",
        "HT55": "HT-55",
        "JM1": "JM-1",
        "K2": "K2 [Human melanoma]",
        "KS-1": "KS-1 [Human glioblastoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "PC-3 [JPC-3]": "PC-3 [Human lung carcinoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "RH-1": "EW-8",
        "SAT": "SAT [Human HNSCC]",
        "TK": "TK [Human B-cell lymphoma]",
    }
    gex = preprocess_df(df=gex, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return gex


def preprocess_methylation_gdsc():
    # read in methylation dataframe
    print("Preprocessing methylation dataframe ...")
    methylation = pd.read_csv("../GDSC/methylation/reprocessed_gdsc_met.csv")
    renamed = {
        "JM1": "JM-1",
        "HT55": "HT-55",
        "MO": "Mo",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "CA-SKI": "Ca Ski",
        "C32": "C32 [Human melanoma]",
        "EOL-1-CELL": "EoL-1",
        "GA-10-CLONE-4": "GA-10 clone 4",
        "HH": "HH [Human lymphoma]",
        "JIYOYEP-2003": "Jiyoye",
        "KS-1": "KS-1 [Human glioblastoma]",
        "LC-2-AD": "LC-2/ad",
        "LNCAP-CLONE-FGC": "LNCaP clone FGC",
        "NO-11": "Onda 11",
        "NO-10": "Onda 10",
        "NTERA-S-CL-D1": "NT2-D1",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "RH-1": "EW-8",
        "HUO-3N1": "HuO-3N1",
        "CAR-1": "CaR-1",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "OMC-1": "OMC-1 [Human cervical carcinoma]",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HEP3B2-1-7": "Hep 3B2.1-7",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "G-292-CLONE-A141B1": "G-292 clone A141B1",
        "HEP_G2": "Hep-G2",
        "LC-1-SQ": "LC-1/sq",
        "RERF-LC-SQ1": "RERF-LC-Sq1",
        "SAT": "SAT [Human HNSCC]",
        "TK": "TK [Human B-cell lymphoma]",
        "PC-3_JPC-3": "PC-3 [Human lung carcinoma]",
    }
    methylation = preprocess_df(df=methylation, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return methylation


def preprocess_mutation_gdsc():
    # read in mutation dataframe
    print("Preprocessing mutation dataframe ...")
    mutation = pd.read_csv("../GDSC/mutation/mutations.csv")
    renamed = {
        "JM1": "JM-1",
        "HT55": "HT-55",
        "K2": "K2 [Human melanoma]",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "C32": "C32 [Human melanoma]",
        "G-292-Clone-A141B1": "G-292 clone A141B1",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HH": "HH [Human lymphoma]",
        "Hep3B2-1-7": "Hep 3B2.1-7",
        "Hs-633T": "Hs 633.T",
        "KS-1": "KS-1 [Human glioblastoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "OMC-1": "OMC-1 [Human cervical carcinoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "SAT": "SAT [Human HNSCC]",
        "TALL-1": "TALL-1 [Human adult T-ALL]",
        "TK": "TK [Human B-cell lymphoma]",
    }
    mutation = preprocess_df(df=mutation, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return mutation


def preprocess_cnv_gdsc():
    # read in copy number variation dataframe
    print("Preprocessing copy number variation dataframe ...")
    cnv = pd.read_csv("../GDSC/cnv/copy_number_variation_gistic.csv", index_col=0)
    renamed = {
        "C32": "C32 [Human melanoma]",
        "G-292-Clone-A141B1": "G-292 clone A141B1",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "Hep3B2-1-7": "Hep 3B2.1-7",
        "HH": "HH [Human lymphoma]",
        "HT55": "HT-55",
        "JM1": "JM-1",
        "K2": "K2 [Human melanoma]",
        "KS-1": "KS-1 [Human glioblastoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "OMC-1": "OMC-1 [Human cervical carcinoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "RH-1": "EW-8",
        "SAT": "SAT [Human HNSCC]",
        "TALL-1": "TALL-1 [Human adult T-ALL]",
        "TK": "TK [Human B-cell lymphoma]",
    }
    cnv = preprocess_df(df=cnv, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return cnv


"""
def preprocess_binarized_drp_gdsc():
    # read in drug response dataframe
    print("Preprocessing drug response dataframe ...")
    drp = pd.read_csv("../GDSC/response/binarized_gdsc.csv")
    drp = drp.drop(0)
    renamed = {
        "JM1": "JM-1",
        "HT55": "HT-55",
        "K2": "K2 [Human melanoma]",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "C32": "C32 [Human melanoma]",
        "SAT": "SAT [Human HNSCC]",
        "TK": "TK [Human B-cell lymphoma]",
        "HH": "HH [Human lymphoma]",
        "G-292 Clone A141B1": "G-292 clone A141B1",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "KS-1": "KS-1 [Human glioblastoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "PC-3 [JPC-3]": "PC-3 [Human lung carcinoma]",
        "OMC-1": "OMC-1 [Human cervical carcinoma]",
        "TALL-1": "TALL-1 [Human adult T-ALL]",
    }
    drp = preprocess_df(df=drp, cl_col_name="Screened Compounds:", renamed_dict=renamed)
    return drp
"""


def preprocess_drp_gdsc(dataset_name):
    # read in drug response dataframe
    print("Preprocessing drug response dataframe ...")
    if dataset_name == "GDSC1":
        path_to_drp = "../GDSC/response/GDSC1_fitted_dose_response_27Oct23.xlsx"
        renamed = {
            "JM1": "JM-1",
            "HT55": "HT-55",
            "K2": "K2 [Human melanoma]",
            "MS-1": "MS-1 [Human lung carcinoma]",
            "C32": "C32 [Human melanoma]",
            "G-292-Clone-A141B1": "G-292 clone A141B1",
            "HARA": "HARA [Human squamous cell lung carcinoma]",
            "HH": "HH [Human lymphoma]",
            "Hep3B2-1-7": "Hep 3B2.1-7",
            "Hs-633T": "Hs 633.T",
            "KS-1": "KS-1 [Human glioblastoma]",
            "ML-1": "ML-1 [Human thyroid carcinoma]",
            "NOS-1": "NOS-1 [Human osteosarcoma]",
            "NTERA-2-cl-D1": "NT2-D1",
            "OMC-1": "OMC-1 [Human cervical carcinoma]",
            "PC-3_[JPC-3]": "PC-3 [Human lung carcinoma]",
            "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
            "RH-1": "EW-8",
            "SAT": "SAT [Human HNSCC]",
            "TALL-1": "TALL-1 [Human adult T-ALL]",
            "TK": "TK [Human B-cell lymphoma]",
        }
    else:
        path_to_drp = "../GDSC/response/GDSC2_fitted_dose_response_27Oct23.xlsx"
        renamed = {
            "HT55": "HT-55",
            "MS-1": "MS-1 [Human lung carcinoma]",
            "C32": "C32 [Human melanoma]",
            "G-292-Clone-A141B1": "G-292 clone A141B1",
            "HARA": "HARA [Human squamous cell lung carcinoma]",
            "HH": "HH [Human lymphoma]",
            "Hep3B2-1-7": "Hep 3B2.1-7",
            "Hs-633T": "Hs 633.T",
            "KS-1": "KS-1 [Human glioblastoma]",
            "K2": "K2 [Human melanoma]",
            "JM1": "JM-1",
            "ML-1": "ML-1 [Human thyroid carcinoma]",
            "NOS-1": "NOS-1 [Human osteosarcoma]",
            "NTERA-2-cl-D1": "NT2-D1",
            "OMC-1": "OMC-1 [Human cervical carcinoma]",
            "PC-3_[JPC-3]": "PC-3 [Human lung carcinoma]",
            "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
            "SAT": "SAT [Human HNSCC]",
            "TALL-1": "TALL-1 [Human adult T-ALL]",
            "TK": "TK [Human B-cell lymphoma]",
        }
    drp = pd.read_excel(path_to_drp)
    drp = drp[["CELL_LINE_NAME", "DRUG_NAME", "LN_IC50", "AUC", "TCGA_DESC", "PUTATIVE_TARGET", "PATHWAY_NAME"]]
    # rename cell line names
    cl_names = drp[["CELL_LINE_NAME"]]
    drp["CELL_LINE_NAME"] = drp["CELL_LINE_NAME"].replace(renamed)
    cl_names = cl_names.drop_duplicates()
    cl_names = preprocess_df(df=cl_names, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return drp, cl_names


def post_process_drp_gdsc(drp, matched_cell_lines, cello_id_to_ac_dict):
    annotation = [matched_cell_lines.get(row['CELL_LINE_NAME']) for index, row in drp.iterrows()]
    drp["cellosaurus_id"] = annotation
    cl_names = [cello_id_to_ac_dict.get(row["cellosaurus_id"]) for index, row in drp.iterrows()]
    drp["cell_line_name"] = cl_names
    drp = drp.drop(columns=["CELL_LINE_NAME"])
    drp = drp.rename(columns={"cell_line_name": "CELL_LINE_NAME"})
    drp = drp.set_index(["cellosaurus_id", "CELL_LINE_NAME"])
    drp = drp.sort_index()
    return drp

def preprocess_sanger_ccle_gex(tpm=True):
    # read in gene expression dataframe
    print("Preprocessing gene expression dataframe ...")
    if tpm:
        gex = pd.read_csv("../SangerCellModelPassports/gene_expression/sanger_tpm_ccle.csv", sep=",")
    else:
        gex = pd.read_csv("../SangerCellModelPassports/gene_expression/sanger_counts_ccle.csv", sep=",")
    return gex


def preprocess_sanger_sanger_gex(tpm=True):
    # read in gene expression dataframe
    print("Preprocessing gene expression dataframe ...")
    if tpm:
        gex = pd.read_csv("../SangerCellModelPassports/gene_expression/sanger_tpm_sanger.csv", sep=",")
    else:
        gex = pd.read_csv("../SangerCellModelPassports/gene_expression/sanger_counts_sanger.csv", sep=",")
    return gex


def preprocess_sanger_proteomcis():
    # read in proteomics dataframe
    print("Preprocessing proteomics dataframe ...")
    prot = pd.read_csv("../SangerCellModelPassports/proteomics/mapped_protein_matrix_maxlfq_diann-normalised.csv", index_col=0)
    id_mapping_df = pd.read_csv('../SangerCellModelPassports/proteomics/uniprot_id_to_gene_name_mapping.tsv', sep='\t')
    prot.index = prot.index.map(id_mapping_df.set_index('from')['to'].to_dict())
    control_runs = prot.columns[prot.columns.str.startswith("Control_HEK293T")].tolist()
    prot = prot.drop(columns=control_runs)
    prot = prot.T
    prot = prot.reset_index(names='cell_line_name')
    # sort by cell line name
    prot = prot.sort_values("cell_line_name")
    # there are 6 replicates per cell line, indicated by .1, .2, …. Remove the .1, .2, …
    prot["cell_line_name"] = prot["cell_line_name"].str.split(".").str[0]
    # kick out all cell_line_names starting with Unnamed
    prot = prot[~prot["cell_line_name"].str.startswith("Unnamed")]
    renamed = {
        "C32": "C32 [Human melanoma]",
        "C3A": "Hep-G2/C3A",
        "G-292-Clone-A141B1": "G-292 clone A141B1",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "Hep3B2-1-7": "Hep 3B2.1-7",
        "HH": "HH [Human lymphoma]",
        "Hs-633T": "Hs 633.T",
        "HT55": "HT-55",
        "JM1": "JM-1",
        "K2": "K2 [Human melanoma]",
        "KS-1": "KS-1 [Human glioblastoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "NTERA-2-cl-D1": "NT2-D1",
        "OMC-1": "OMC-1 [Human cervical carcinoma]",
        "PC-3_[JPC-3]": "PC-3 [Human lung carcinoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "RH-1": "EW-8",
        "SAT": "SAT [Human HNSCC]",
        "TALL-1": "TALL-1 [Human adult T-ALL]",
        "TK": "TK [Human B-cell lymphoma]",
        "UWB1": "UWB1.289",

    }
    prot = preprocess_df(df=prot, cl_col_name="cell_line_name", renamed_dict=renamed)
    return prot


def preprocess_ccle_gex():
    # read in gene expression dataframe
    print("Preprocessing gene expression dataframe ...")
    gex = pd.read_csv("../CCLE/gene_expression/gene_expression.csv")
    renamed = {
        "C32": "C32 [Human melanoma]",
        "RCM-1": "RCM-1 [Human rectal adenocarcinoma]",
        "PC-3 [JPC-3]": "PC-3 [Human lung carcinoma]",
        "RH-1": "EW-8",
        "KS-1": "KS-1 [Human glioblastoma]",
        "G-292 Clone A141B1": "G-292 clone A141B1",
        "MS-1": "MS-1 [Human lung carcinoma]",
        "ML-1": "ML-1 [Human thyroid carcinoma]",
        "SAT": "SAT [Human HNSCC]",
        "HH": "HH [Human lymphoma]",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "TK": "TK [Human B-cell lymphoma]",
        "K2": "K2 [Human melanoma]",
        "NOS-1": "NOS-1 [Human osteosarcoma]",
        "HT55": "HT-55",
        "JM1": "JM-1",
    }
    gex = preprocess_df(df=gex, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return gex


def preprocess_drp_CCLE():
    print("Preprocessing CCLE response dataframe ...")
    drp = pd.read_csv("../CCLE/response/response_CCLE.csv")
    drp = drp.rename(columns={"CELL_LINE_NAME": "cell_line_name"})
    # get long format into wide format: rows should be cell line names (CELL_LINE_NAME column), columns should be drug names (DRUG_NAME column), values are in LN_IC50 column
    drp = drp.pivot(index="cell_line_name", columns="DRUG_NAME", values="LN_IC50")
    drp = drp.reset_index()
    renamed = {
        "C32": "C32 [Human melanoma]",
        "C3A": "Hep-G2/C3A",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HH": "HH [Human lymphoma]",
        "JM1": "JM-1",
    }
    drp = preprocess_df(df=drp, cl_col_name="cell_line_name", renamed_dict=renamed)
    return drp


def preprocess_methylation_ccle():
    print("Preprocessing CCLE methylation dataframe ...")
    met = pd.read_csv("../CCLE/methylation/methylation_pre_mapping.csv")
    # cell line name has tissue appended to it, separated by "_". Remove tissue information.
    met["CELL_LINE"] = met["CELL_LINE"].str.split("_").str[0]
    renamed = {
        "BC3C": "BC-3C",
        "C32": "C32 [Human melanoma]",
        "CJM": "CJM [Human melanoma]",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HH": "HH [Human lymphoma]",
        "JK1": "JK-1",
        "JM1": "JM-1",
        "KG1": "KG-1",
        "ML1": "ML-1 [Human thyroid carcinoma]",
        "NCIH2087": "NCI-H2087",
        "NCIH2444": "NCI-H2444",
        "PC3": "PC-3",
        "PK1": "PK-1",
        "RCM1": "RCM-1 [Human rectal adenocarcinoma]",
        "SH4": "SH-4",
        "TE1": "TE-1",
        "TE8": "TE-8",
        "TE15": "TE-15",
    }
    met = preprocess_df(df=met, cl_col_name="CELL_LINE", renamed_dict=renamed)
    return met


def preprocess_proteomics_ccle():
    print("Preprocessing CCLE proteomics dataframe ...")
    prot = pd.read_csv("../CCLE/proteomics/protein_expression_table.tsv", sep="\t")
    prot = prot.set_index("Gene names")
    prot = prot.T
    prot = prot.reset_index(names='cell_line_name')
    renamed = {
        "C32": "C32 [Human melanoma]",
        "CAL-120.1": "CAL-120",
        "CNI-H1568": "NCI-H1568",
        "HCC 1806": "HCC1806",
        "HCC 1954": "HCC1954",
        "HCT-15.1": "HCT 15",
        "HT55": "HT-55",
        "JM1": "JM-1",
        "SW948.1": "SW948",
    }
    prot = preprocess_df(df=prot, cl_col_name="cell_line_name", renamed_dict=renamed)
    return prot


def preprocess_reprocessed_ccle_gistic():
    print("Preprocessing reprocessed CCLE GISTIC dataframe ...")
    cnv = pd.read_csv("../CCLE/cnv/reprocessed_cnv.csv")
    renamed = {
        "BC3C": "BC-3C",
        "C32": "C32 [Human melanoma]",
        "C3A": "Hep-G2/C3A",
        "CJM": "CJM [Human melanoma]",
        "COLO775": "COLO 775",
        "HARA": "HARA [Human squamous cell lung carcinoma]",
        "HH": "HH [Human lymphoma]",
        "HK2": "HK-2 [Human kidney]",
        "HT29": "HT-29",
        "HT55": "HT-55",
        "HTK": "HTK-",
        "JM1": "JM-1",
        "JK1": "JK-1",
        "KG1": "KG-1",
        "MEC1": "MEC-1",
        "ML1": "ML-1 [Human thyroid carcinoma]",
        "NCIH2087": "NCI-H2087",
        "NCIH2444": "NCI-H2444",
        "PC3": "PC-3",
        "PK1": "PK-1",
        "RCM1": "RCM-1 [Human rectal adenocarcinoma]",
        "SH4": "SH-4",
        "TE1": "TE-1",
        "TE8": "TE-8",
        "TE15": "TE-15",
        "U178": "U-178MG",
    }
    cnv = preprocess_df(df=cnv, cl_col_name="CELL_LINE_NAME", renamed_dict=renamed)
    return cnv


if __name__ == "__main__":
    cellosaurus = pd.read_csv("../mapping/cellosaurus_01_2024.csv")
    # replace all NaN values with empty strings
    cellosaurus = cellosaurus.fillna("")
    # create cellosaurus dictionary
    cellosaurus_ac_dict, cellosaurus_sy_dict, species_dict, cello_ac_to_id_dict = create_cl_dict(cellosaurus)

    # GDSC
    # map gene expression cell line names to cellosaurus IDs
    gex = preprocess_gex_gdsc()
    map_to_cellosaurus(
        gex,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../GDSC/gene_expression/gene_expression_cellosaurus.csv",
    )
    
    # map methylation cell line names to cellosaurus IDs
    met = preprocess_methylation_gdsc()
    map_to_cellosaurus(
        met,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../GDSC/methylation/methylation_cellosaurus.csv",
    )

    # map copy number variation cell line names to cellosaurus IDs
    cnv = preprocess_cnv_gdsc()
    map_to_cellosaurus(
        cnv,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../GDSC/cnv/copy_number_variation_gistic_cellosaurus.csv",
    )
    
    # map drug response cell line names to cellosaurus IDs
    drp_gdsc_1, cl_names = preprocess_drp_gdsc(dataset_name="GDSC1")
    map_to_cellosaurus(
        cl_names,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../mapping/cl_mapping_gdsc1.csv",
    )
    drp_gdsc_1 = post_process_drp_gdsc(drp_gdsc_1, matched_cell_lines, cello_ac_to_id_dict)
    drp_gdsc_1.to_csv("../GDSC/response/response_GDSC1_cellosaurus.csv", index=True)

    # map drug response cell line names to cellosaurus IDs
    drp_gdsc_2, cl_names = preprocess_drp_gdsc(dataset_name="GDSC2")
    map_to_cellosaurus(
        cl_names,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../mapping/cl_mapping_gdsc2.csv",
    )
    drp_gdsc_2 = post_process_drp_gdsc(drp_gdsc_2, matched_cell_lines, cello_ac_to_id_dict)
    drp_gdsc_2.to_csv("../GDSC/response/response_GDSC2_cellosaurus.csv", index=True)
    
    #cellosaurus_sid_dict, species_dict = create_cl_dict_cell_passp(cellosaurus)
    # CCLE data
    # response
    ccle_drp = preprocess_drp_CCLE()
    map_to_cellosaurus(
        ccle_drp,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../CCLE/response/response_long_CCLE_cellosaurus.csv",
    )
    ccle_drp = ccle_drp.melt(id_vars=["cellosaurus_id", "cell_line_name"], var_name="DRUG_NAME", value_name="LN_IC50")
    # rename columns: cell_line_name -> CELL_LINE_NAME
    ccle_drp = ccle_drp.rename(columns={"cell_line_name": "CELL_LINE_NAME"})
    ccle_drp = ccle_drp.set_index("cellosaurus_id")
    ccle_drp.to_csv("../CCLE/response/response_CCLE_cellosaurus.csv", index=True)

    # gene expression
    ccle = preprocess_ccle_gex()
    map_to_cellosaurus(
        ccle,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../CCLE/gene_expression/gene_expression_cellosaurus.csv",
    )
    # methylation
    ccle_met = preprocess_methylation_ccle()
    map_to_cellosaurus(
        ccle_met,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../CCLE/methylation/methylation_cellosaurus.csv",
    )

    # reprocessed GISTIC data
    ccle_cnv = preprocess_reprocessed_ccle_gistic()
    map_to_cellosaurus(
        ccle_cnv,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../CCLE/cnv/cnv_cellosaurus.csv",
    )
    
    # proteomics
    ccle_prot = preprocess_proteomics_ccle()
    map_to_cellosaurus(
        ccle_prot,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../CCLE/proteomics/proteomics_cellosaurus.csv",
    )

    sanger_prot = preprocess_sanger_proteomcis()
    map_to_cellosaurus(
        sanger_prot,
        cellosaurus_ac_dict,
        cellosaurus_sy_dict,
        species_dict,
        cello_ac_to_id_dict,
        "../SangerCellModelPassports/proteomics/proteomics_cellosaurus.csv",
    )

    # export matched cell lines to csv file
    matched_cell_lines_df = pd.DataFrame.from_dict(matched_cell_lines, orient="index", columns=["cellosaurus_id"])
    matched_cell_lines_df.to_csv("../mapping/matched_cell_lines.csv")

    # export unmatched cell lines to csv file
    unmatched_cell_lines_df = pd.DataFrame(list(no_match), columns=["cell_line_name"])
    unmatched_cell_lines_df.to_csv("../mapping/unmatched_cell_lines.csv")

    # export renamed cell lines to csv file
    renamed_cell_lines_df = pd.DataFrame.from_dict(renamed_cell_lines, orient="index", columns=["cell_line_name"])
    renamed_cell_lines_df.to_csv("../mapping/renamed_cell_lines.csv")

