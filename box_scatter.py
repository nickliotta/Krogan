import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

data_set = pd.read_csv("data/results-mss-normalized.csv")

custom_palette = {
    "KO.KO.1" : "#df6464",
    "KO.KO.2" : "#df6464",
    "KO.KO.3" : "#df6464",
    "KO.KO.4" : "#df6464",
    "WT.WT.1" : "#3d85c6",
    "WT.WT.2" : "#3d85c6",
    "WT.WT.3" : "#3d85c6",
    "WT.WT.4" : "#3d85c6"
}

sns.set_theme(style="whitegrid")
flierprops = dict(marker="o", markerfacecolor="black", markersize=3,  markeredgecolor="black")

sns.boxplot(
    data=data_set, 
    x=data_set["SUBJECT_ORIGINAL"], 
    y=np.log2(data_set["INTENSITY"]),
    palette=custom_palette,
    flierprops=flierprops
)

cp = {
    "KO.KO.1" : "#cc0000",
    "KO.KO.2" : "#cc0000",
    "KO.KO.3" : "#cc0000",
    "KO.KO.4" : "#cc0000",
    "WT.WT.1" : "#0b5394",
    "WT.WT.2" : "#0b5394",
    "WT.WT.3" : "#0b5394",
    "WT.WT.4" : "#0b5394"
}

sns.stripplot(
    data=data_set, 
    x=data_set["SUBJECT_ORIGINAL"], 
    y=np.log2(data_set["INTENSITY"]),
    size=1,
    marker="s", 
    palette=cp,
    alpha=0.2
)

plt.title("Protein Intensity in BioReplicates (Excluding contaminants)")
plt.suptitle("")
plt.ylabel("Intensity")
plt.xticks(rotation=90)
plt.xlabel("BioReplicate")
plt.show()