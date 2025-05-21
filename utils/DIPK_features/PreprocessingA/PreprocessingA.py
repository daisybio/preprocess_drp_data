# code adapted from https://github.com/reinej03/DIPK_Bachelor_Thesis/tree/main and https://github.com/user15632/DIPK
import numpy as np
import pandas as pd

gene_add_num = 256

#gene_list_sel.py -------------------------------------------------------------------------------------------------------------------------
"""
f = open('../Data/key.genes.txt', encoding='gbk')
gene_list_1 = []
for each_row in f:
    gene_list_1.append(each_row.strip())

dataset = pd.read_csv('../Data/Cell_line_RMA_proc_basalExp.txt', header=None, sep='\t', low_memory=False)
gene_list_2 = [str(_) for _ in list(dataset.iloc[1:, 0])]

genes = sorted(list(set(gene_list_1) & set(gene_list_2)))

for each in genes:
    with open('../Data/gene_list_sel.txt', 'a') as file0:
        print(each, file=file0) 
"""
#RMA_dict.py -------------------------------------------------------------------------------------------------------------------------
print("Running: RMA")

f = open('../Data/gene_list_sel.txt', encoding='gbk')
gene_list = []
for each_row in f:
    gene_list.append(each_row.strip())

dataset = pd.read_csv('../Data/Cell_line_RMA_proc_basalExp.txt', header=None, sep='\t', low_memory=False)
Cell = ['.'.join(str(_).split('.')[1:]) for _ in list(dataset.iloc[0, 2:])]
Gene = [str(_) for _ in list(dataset.iloc[1:, 0])]
Gene_index = [Gene.index(_) for _ in gene_list]

dataset = pd.read_csv('../Data/Cell_line_RMA_proc_basalExp.txt', header=0, sep='\t')
RMA_df = dataset.iloc[:, 2:]
RMA_df = RMA_df.iloc[Gene_index, :]
RMA_df = RMA_df.T
RMA_df.index=Cell
RMA_df.columns = range(RMA_df.shape[1])

RMA_df.to_csv('../Data/RMA.csv', sep='\t')

#BIONIC.py -------------------------------------------------------------------------------------------------------------------------
print("Running: BIONC")

f = open('../Data/gene_list_sel.txt', encoding='gbk')
gene_list = []
for each_row in f:
    gene_list.append(each_row.strip())

dataset = pd.read_csv('../Data/human_ppi_features.tsv', header=0, sep='\t', index_col=0)
index = [i for i, elem in enumerate(gene_list) if elem in dataset.index] #this line was taken from ChatGPT! [CG1]
BIONIC_df=dataset[dataset.index.isin(gene_list)]
BIONIC_df.index=index

BIONIC_df.to_csv('../Data/BIONIC.csv', sep='\t')

#Data.py -------------------------------------------------------------------------------------------------------------------------
indices = np.argsort(-(RMA_df.values))

features=[]
for i in range(RMA_df.shape[0]):
    k = 0
    current_feature = np.zeros(BIONIC_df.shape[1])
    for j in range(gene_add_num):
        if indices[i, j] in BIONIC_df.index:
            k += 1
            current_feature += BIONIC_df.loc[indices[i,j]]
    current_feature = current_feature / k
    features.append(current_feature)

BNF_df=pd.DataFrame(features, index=RMA_df.index)

#Harmonize Data --------------------------------------------------------------------------------------------------------------------
print("Harmonizing Data")
cell_line_details = pd.read_csv('../Data/Cell_Lines_Details.csv', header=0, sep=';', index_col=1)

cell_line_details.index = cell_line_details.index.astype(str)
BNF_df = BNF_df[BNF_df.index.isin(cell_line_details.index)]
BNF_df.index = cell_line_details.loc[BNF_df.index, 'Line']

BNF_df.to_csv('../Data/BNF.csv', sep='\t')

print(BNF_df)





