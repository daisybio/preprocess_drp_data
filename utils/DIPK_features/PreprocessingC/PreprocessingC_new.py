# code adapted from https://github.com/reinej03/DIPK_Bachelor_Thesis/tree/main and https://github.com/user15632/DIPK
import pandas as pd
import joblib
from rdkit import Chem
import numpy as np
import torch
from torch_geometric.data import Data
import os

from PreprocessingC_loader import mol_to_graph_data_obj_complex

from PreprocessingC_Model_GNN import *
from PreprocessingC_util import *


#Smiles ------------------------------------------------------------------------------------------------------------------------------
import pandas as pd
import os
os
all_smiles = pd.read_csv('drugs/all_smiles_final.csv')
all_smiles.dropna(subset=['canonical_smiles'], inplace=True)
SMILES_dict = {pbid: sm for pbid, sm in zip(all_smiles['pubchem_id'], all_smiles['canonical_smiles'])}
#Graphs ----------------------------------------------------------------------------------------------------------------------------
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
gnn.load_state_dict(torch.load('MolGNet.pt'))
gnn = gnn.to(DEVICE)
gnn.eval()
with torch.no_grad():
    for each in GRAPH_dict:
        print(each)
        graph = GRAPH_dict[each]
        graph = Self_loop(graph)
        graph = Add_seg_id(graph)
        graph = graph.to(DEVICE)
        MolGNet_dict[each] = gnn(graph).cpu()
joblib.dump(MolGNet_dict, '../Data/MolGNet_dict.pkl')

# CSVs -----------------------------------------------------------------------------------------------------------------------------
MolGNet_dict = joblib.load('../Data/MolGNet_dict.pkl')
for each in MolGNet_dict.keys():
    MolGNet_df = pd.DataFrame.from_dict(MolGNet_dict[each])
    MolGNet_df.to_csv(f"../Data/Drugs_new/MolGNet_{each}.csv", sep='\t')
    
if False:
    GRAPH_dict = joblib.load('../Data/GRAPH_dict.pkl')
    for each in GRAPH_dict:
        graph = GRAPH_dict[each]
        GRAPH_df = pd.DataFrame.from_dict(graph.edge_index)
        GRAPH_df.to_csv(f"../Data/Drugs_new/{each}/Edge_Index_{each}.csv", sep='\t')  
    for each in GRAPH_dict:
        graph = GRAPH_dict[each]
        GRAPH_df = pd.DataFrame.from_dict(graph.edge_attr)
        GRAPH_df.to_csv(f"../Data/Drugs_new/{each}/Edge_Attr_{each}.csv", sep='\t')

