from bioinfokit import visuz
import pandas as pd

data_set = pd.read_csv("data/results.csv")
visuz.GeneExpression.volcano(
    df=data_set, 
    lfc="log2FC", 
    pv="pvalue", 
    show=True,
    color=("#ff0000", "grey", "red"),
    sign_line=True,
    axxlabel="log2FC",
    axylabel="-log10(pvalue)"
)