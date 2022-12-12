from matplotlib.patches import Rectangle
from matplotlib import axes, colors
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

_, axes = plt.subplots()
axes.yaxis.tick_right()

raw_df = pd.read_csv("data/heatmap/test4.csv")
df = raw_df.drop("pval", axis=1)
df = df.pivot("label", "term", "nes")
df = df.transpose()

cmap = colors.LinearSegmentedColormap.from_list("n", [
    "#001eff",
    "#ffffff", # middle
    "#cc0000"])

ax = sns.heatmap(
        df, 
        vmin=-4, 
        vmax=4,
        xticklabels=True, 
        yticklabels=True,
        fmt=".2f", 
        cmap=cmap, 
        linewidths=0.3, 
        linecolor="gray",
        cbar_kws={"shrink": .40, "pad": 0.03, "location" : "left", "label": "NES", "extend" : "max"})

for index, row in raw_df.iterrows():
    pval = row[3]
    if pval <= 0.05:
        row_ = [i for i, x in enumerate(df.index == row["term"]) if x][0]
        col = [i for i, x in enumerate(df.columns == row["label"]) if x][0]

        ax.add_patch(Rectangle((col, row_), 1, 1, fill=False, edgecolor="black", lw=0.8))

ax.tick_params(length=0.5)
ax.set_aspect("equal")
ax.set_ylabel("")    
ax.set_xlabel("")

plt.title(label="Gene Ontology ( Top 30 )", loc="left", pad=10)
plt.subplots_adjust(bottom=0.23)
plt.xticks(rotation=90, fontsize=8)
plt.yticks(fontsize=8)
plt.show()


