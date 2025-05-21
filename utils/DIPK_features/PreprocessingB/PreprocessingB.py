# code adapted from https://github.com/reinej03/DIPK_Bachelor_Thesis/tree/main and https://github.com/user15632/DIPK

import pandas as pd 
import torch
import joblib

#feature extraction with DAE ------------------------------------------------------------------------------------------------------
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

encoder, _ = joblib.load('PreTrain.pkl')
encoder = encoder.to(DEVICE)

RMA_df = pd.read_csv('../Data/RMA.csv', sep='\t', index_col=0)
RMA_tensor = torch.tensor(RMA_df.values, dtype=torch.float32).to(DEVICE)

encoder.eval()
with torch.no_grad():
    GEF = encoder(RMA_tensor)

GEF_df = pd.DataFrame(GEF, index=RMA_df.index)


#Harmonize Data --------------------------------------------------------------------------------------------------------------------
cell_line_details = pd.read_csv('../Data/Cell_Lines_Details.csv', header=0, sep=';', index_col=1)

cell_line_details.index = cell_line_details.index.astype(str)
GEF_df.index = GEF_df.index.astype(str).str.replace(r'\.0', '', regex=True) #this line was taken from ChatGPT! [CG2]
GEF_df = GEF_df[GEF_df.index.isin(cell_line_details.index)]
GEF_df.index = cell_line_details.loc[GEF_df.index, 'Line']

GEF_df.to_csv('../Data/GEF.csv')




