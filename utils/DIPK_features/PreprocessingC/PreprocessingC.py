# code adapted from https://github.com/reinej03/DIPK_Bachelor_Thesis/tree/main and https://github.com/user15632/DIPK
import pandas as pd
import pubchempy as pcp
import joblib
from rdkit import Chem
import numpy as np
import torch
from torch_geometric.data import Data
import os

from PreprocessingC_loader import mol_to_graph_data_obj_complex

from PreprocessingC_Model_GNN import *
from PreprocessingC_util import *

#Drugnames -------------------------------------------------------------------------------------------------------------------------
train_set = pd.read_csv('../Data/fold_0_train_cell_drug_pair.csv', header=None)
Ntrain = list(train_set.iloc[:, 1])
val_set = pd.read_csv('../Data/fold_0_validation_cell_drug_pair.csv', header=None)
Nval = list(val_set.iloc[:, 1])
test_set = pd.read_csv('../Data/fold_0_test_cell_drug_pair.csv', header=None)
Ntest = list(test_set.iloc[:, 1])

drug_names = Ntrain + Nval + Ntest 

drug_names = set(drug_names)
drug_names = list(drug_names)

#Smiles ------------------------------------------------------------------------------------------------------------------------------
pubchem_cid = pd.read_csv('../Data/pubchem_cid.csv', header=0)
pubchem_cid = pubchem_cid[~(pubchem_cid.iloc[:, 5].isin(['none', 'several']))]
pubchem_cid = pubchem_cid[pubchem_cid.iloc[:, 5].notna()]

name = list(pubchem_cid.iloc[:, 1])
cid = list(pubchem_cid.iloc[:, 5])
pubchem_dict = dict(zip(name, cid))

SMILES_dict = dict()
SMILES_dict['Lestauritinib'] = 'CC12C(CC(O1)N3C4=CC=CC=C4C5=C6C(=C7C8=CC=CC=C8N2C7=C53)CNC6=O)(CO)O'

SMILES_not_found = []
i = 0
for each in drug_names:
    i += 1
    try:
        if each in pubchem_dict:
            _ = pcp.Compound.from_cid(pubchem_dict[each])
            SMILES_dict[each] = _.isomeric_smiles
        else:
            _ = pcp.get_compounds(each, 'name')
            SMILES_dict[each] = _[0].isomeric_smiles
    except:
        SMILES_not_found.append(each)
        print(each)
    print(i, end=' ')
print()
print("Smiles not found:", SMILES_not_found)  # ['KIN001-236', 'BX796']
print("Length of Smiles not found:", len(SMILES_not_found))
joblib.dump(SMILES_dict, '../Data/SMILES_dict.pkl')
#Graphs ----------------------------------------------------------------------------------------------------------------------------
SMILES_dict = joblib.load('../Data/SMILES_dict.pkl')
GRAPH_dict = dict()
for each in SMILES_dict:
    GRAPH_dict[each] = mol_to_graph_data_obj_complex(Chem.MolFromSmiles(SMILES_dict[each]))
    # Data(x=[num_nodes, 8], edge_index=[2, num_edges], edge_attr=[num_edges, 5])
joblib.dump(GRAPH_dict, '../Data/GRAPH_dict.pkl')

#MolGNet ---------------------------------------------------------------------------------------------------------------------------
GRAPH_dict = joblib.load('../Data/GRAPH_dict.pkl')
MolGNet_dict = dict()

Self_loop = Self_loop()
Add_seg_id = Add_seg_id()

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
gnn = MolGNet(num_layer=5, emb_dim=768, heads=12, num_message_passing=3, drop_ratio=0)
gnn.load_state_dict(torch.load('PreprocessingC_MolGNet.pt'))
gnn = gnn.to(DEVICE)
gnn.eval()
with torch.no_grad():
    for each in GRAPH_dict:
        graph = GRAPH_dict[each]
        graph = Self_loop(graph)
        graph = Add_seg_id(graph)
        graph = graph.to(DEVICE)
        MolGNet_dict[each] = gnn(graph).cpu()
joblib.dump(MolGNet_dict, '../Data/MolGNet_dict.pkl')

# CSVs -----------------------------------------------------------------------------------------------------------------------------
MolGNet_dict = joblib.load('../Data/MolGNet_dict.pkl')
for each in MolGNet_dict.keys():
    os.makedirs(f"../Data/Drugs/{each}", exist_ok=True) 
    MolGNet_df = pd.DataFrame.from_dict(MolGNet_dict[each])
    MolGNet_df.to_csv(f"../Data/Drugs/{each}/MolGNet_{each}.csv", sep='\t')
    
GRAPH_dict = joblib.load('../Data/GRAPH_dict.pkl')
for each in GRAPH_dict:
    graph = GRAPH_dict[each]
    GRAPH_df = pd.DataFrame.from_dict(graph.edge_index)
    GRAPH_df.to_csv(f"../Data/Drugs/{each}/Edge_Index_{each}.csv", sep='\t')  
for each in GRAPH_dict:
    graph = GRAPH_dict[each]
    GRAPH_df = pd.DataFrame.from_dict(graph.edge_attr)
    GRAPH_df.to_csv(f"../Data/Drugs/{each}/Edge_Attr_{each}.csv", sep='\t')

