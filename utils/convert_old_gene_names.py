import pandas as pd

gene_list = pd.read_csv("../gene_lists/landmark_genes.csv")

new_gene_names = {
    "AARS": "AARS1",
    "EPRS": "EPRS1",
    "FAM57A": "TLCD3A",
    "FAM69A": "DIPK1A",
    "H2AFV": "H2AZ2",
    "HIST1H2BK": "H2BC12",
    "HIST2H2BE": "H2BC21",
    "KIAA0100": "BLTP2",
    "KIAA0355": "GARRE1",
    "NARFL": "CIAO3",
    "PAPD7": "TENT4A",
    "SKIV2L": "SKIC2",
    "TSTA3": "GFUS",
    "WDR61": "SKIC8",
    "WRB": "GET1",
}
"""
gene_list["Symbol"] = gene_list["Symbol"].replace(new_gene_names)
gene_list.to_csv("../gene_lists/landmark_genes.csv", index=False)
"""
omics = pd.read_csv("copy_number_variation_gistic.csv")
# replace old gene names with new gene names in omics columns
omics = omics.rename(columns=new_gene_names)
genes_in_list = set(gene_list["Symbol"])
genes_in_features = set(omics.columns)
missing_genes = list(genes_in_list - genes_in_features)
missing_genes.sort()
print(missing_genes)