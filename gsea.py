from gseapy.plot import gseaplot, heatmap
import pandas as pd
import gseapy as gp
import numpy as np

names = pd.read_csv("./data/names.txt", delimiter="\t")
def protein_to_gene_name(protein: str):
    return names.loc[names["Entry"] == protein]["Gene_names"].to_string(index=False)

def data_set_to_rnk(name: str):
    df = pd.read_csv(f"./data/{name}.csv")
    cleaned = []

    for _, row in df.iterrows():
        proteins = row["Protein"].split(";")
        for protein in proteins:
            try:
                protein = protein[0:protein.index("_")]
            except ValueError:
                pass

            gene = protein_to_gene_name(protein)

            cleaned.append([
                protein,
                gene,
                np.abs(row["log2FC"]),
                row["log2FC"]    
            ])

    df = pd.DataFrame(cleaned, columns=["Protein", "Gene_name", "abs(log2FC)", "log2FC"])
    df = df.sort_values(by="abs(log2FC)", ascending=False)
    df = df.drop_duplicates(subset="Gene_name", keep="first")

    df.drop(["Protein", "abs(log2FC)"], inplace=True, axis=1)
    df.sort_values(by="log2FC", ascending=False)
    df.to_csv(f"{name}.rnk", sep="\t", header=None, index=False)

rnk = pd.read_csv("./gsea/ubiquitination/ubiquitination.rnk", header=None, delimiter="\t")
pre_res = gp.prerank(rnk=rnk,
                     gene_sets="./data/go.symbols.gmt",
                     processes=4,
                     permutation_num=100000)

pre_res.res2d.to_csv("./gsea/ubiquitination/ubiquitination_go.csv")
