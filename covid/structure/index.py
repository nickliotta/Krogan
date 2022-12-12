
# libraries
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import matplotlib
import pandas as pd 
import json
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import subprocess

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

font_properties = fm.FontProperties(fname="/Users/nick/miniconda/pkgs/matplotlib-base-3.2.2-py38h1300a51_1/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/Ubuntu-Light.ttf")

with open("data.json") as file:
    data = json.load(file)

kwargs = {
    "edgecolor" : "black",
    "lw" : 1,
    "height" : .5,
    "zorder" : 2
}

variants = {
    "omicron" : 0,
    "gamma" : 4,
    "delta" : 8,
    "beta" : 12,
    "alpha" : 16
}

df = pd.read_csv("data.csv")
formulas = {
    "start" : lambda position: -1 + (position * 0.01),
    "ORF1a" : lambda position: 1.65 + ((position*3) * 0.01),
    "ORF1b" : lambda position: 133.8 + ((position*3) * 0.01),
    "S" : lambda position: 214.72 + ((position*3) * 0.01),
    "ORF3a" : lambda position: 253.02 + ((position*3) * 0.01),   
    "E" : lambda position: 261.54 + ((position*3) * 0.01),        
    "M" : lambda position: 264.32 + ((position*3) * 0.01),       
    "ORF7a" : lambda position: 273.03 + ((position*3) * 0.01),       
    "ORF7b" : lambda position: 276.69 + ((position*3) * 0.01),    
    "ORF8" : lambda position: 278.07 + ((position*3) * 0.01),      
    "N" : lambda position: 281.87 + ((position*3) * 0.01),        
    "ORF9b" : lambda position: 281.97 + ((position*3) * 0.01),      
    "ORF10" : lambda position: 294.7 + ((position*3) * 0.01),       
}


fig, ax = plt.subplots(figsize=(16, 10), dpi=80)

ax.vlines(x=1, ymin=0, ymax=5, color="white", alpha=0.7, linewidth=2)
ax.scatter(x=1, y=18, s=75, color="white", alpha=0.7)

ax.set_xticks([])
ax.set_yticks([])

for variant in variants.items():
    # 0.01 per nucleotides
    # Model done in nucleic acids
    # Positions of nonsynonymous mutations are multiplied by 3x to account for nucleic acid position (3 in 1)
    # Positions of synonymous mutations are not changed, and start from the beginning of the genome

    # n/a
    plt.gca().add_patch(patches.Rectangle((-1, variant[1]), width=2.65  * 1.4, facecolor="#e0dcdb", **kwargs))

    # ORF1a
    ORF1a = patches.Rectangle((1.65 * 1.4, variant[1]), width=132.15 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a)

    ORF1a_rx, ORF1a_ry = ORF1a.get_xy()
    ORF1a_cx = ORF1a_rx + ORF1a.get_width()/2.0
    ORF1a_cy = ORF1a_ry + ORF1a.get_height()/2.0

    ORF1a_nsp1 = patches.Rectangle((1.65 * 1.4, variant[1]), width=5.4 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp1)

    ORF1a_nsp1_rx, ORF1a_nsp1_ry = ORF1a_nsp1.get_xy()
    ORF1a_nsp1_cx = ORF1a_nsp1_rx + ORF1a_nsp1.get_width()/2.0
    ORF1a_nsp1_cy = ORF1a_nsp1_ry + ORF1a_nsp1.get_height()/2.0

    ORF1a_nsp2 = patches.Rectangle((7.05 * 1.4, variant[1]), width=19.14 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp2)

    ORF1a_nsp2_rx, ORF1a_nsp2_ry = ORF1a_nsp2.get_xy()
    ORF1a_nsp2_cx = ORF1a_nsp2_rx + ORF1a_nsp2.get_width()/2.0
    ORF1a_nsp2_cy = ORF1a_nsp2_ry + ORF1a_nsp2.get_height()/2.0

    ORF1a_nsp3 = patches.Rectangle((26.19 * 1.4, variant[1]), width=58.35 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp3)

    ORF1a_nsp3_rx, ORF1a_nsp3_ry = ORF1a_nsp3.get_xy()
    ORF1a_nsp3_cx = ORF1a_nsp3_rx + ORF1a_nsp3.get_width()/2.0
    ORF1a_nsp3_cy = ORF1a_nsp3_ry + ORF1a_nsp3.get_height()/2.0

    ORF1a_nsp4 = patches.Rectangle((84.54 * 1.4, variant[1]), width=15 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp4)

    ORF1a_nsp4_rx, ORF1a_nsp4_ry = ORF1a_nsp4.get_xy()
    ORF1a_nsp4_cx = ORF1a_nsp4_rx + ORF1a_nsp4.get_width()/2.0
    ORF1a_nsp4_cy = ORF1a_nsp4_ry + ORF1a_nsp4.get_height()/2.0

    ORF1a_nsp5 = patches.Rectangle((99.54 * 1.4, variant[1]), width=9.18 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp5)

    ORF1a_nsp5_rx, ORF1a_nsp5_ry = ORF1a_nsp5.get_xy()
    ORF1a_nsp5_cx = ORF1a_nsp5_rx + ORF1a_nsp5.get_width()/2.0
    ORF1a_nsp5_cy = ORF1a_nsp5_ry + ORF1a_nsp5.get_height()/2.0

    ORF1a_nsp6 = patches.Rectangle((108.72 * 1.4, variant[1]), width=8.7 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp6)

    ORF1a_nsp6_rx, ORF1a_nsp6_ry = ORF1a_nsp6.get_xy()
    ORF1a_nsp6_cx = ORF1a_nsp6_rx + ORF1a_nsp6.get_width()/2.0
    ORF1a_nsp6_cy = ORF1a_nsp6_ry + ORF1a_nsp6.get_height()/2.0

    ORF1a_nsp7 = patches.Rectangle((117.42 * 1.4, variant[1]), width=2.49 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp7)

    ORF1a_nsp7_rx, ORF1a_nsp7_ry = ORF1a_nsp7.get_xy()
    ORF1a_nsp7_cx = ORF1a_nsp7_rx + ORF1a_nsp7.get_width()/2.0
    ORF1a_nsp7_cy = ORF1a_nsp7_ry + ORF1a_nsp7.get_height()/2.0

    ORF1a_nsp8 = patches.Rectangle((119.91 * 1.4, variant[1]), width=5.94 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp8)

    ORF1a_nsp8_rx, ORF1a_nsp8_ry = ORF1a_nsp8.get_xy()
    ORF1a_nsp8_cx = ORF1a_nsp8_rx + ORF1a_nsp8.get_width()/2.0
    ORF1a_nsp8_cy = ORF1a_nsp8_ry + ORF1a_nsp8.get_height()/2.0

    ORF1a_nsp9 = patches.Rectangle((125.85 * 1.4, variant[1]), width=3.39 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp9)

    ORF1a_nsp9_rx, ORF1a_nsp9_ry = ORF1a_nsp9.get_xy()
    ORF1a_nsp9_cx = ORF1a_nsp9_rx + ORF1a_nsp9.get_width()/2.0
    ORF1a_nsp9_cy = ORF1a_nsp9_ry + ORF1a_nsp9.get_height()/2.0

    ORF1a_nsp10 = patches.Rectangle((129.24 * 1.4, variant[1]), width=4.17 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp10)

    ORF1a_nsp10_rx, ORF1a_nsp10_ry = ORF1a_nsp10.get_xy()
    ORF1a_nsp10_cx = ORF1a_nsp10_rx + ORF1a_nsp10.get_width()/2.0
    ORF1a_nsp10_cy = ORF1a_nsp10_ry + ORF1a_nsp10.get_height()/2.0

    ORF1a_nsp11 = patches.Rectangle((133.41 * 1.4, variant[1]), width=0.39 * 1.4, facecolor="#a7c5df", **kwargs)
    plt.gca().add_patch(ORF1a_nsp11)

    ORF1a_nsp11_rx, ORF1a_nsp11_ry = ORF1a_nsp11.get_xy()
    ORF1a_nsp11_cx = ORF1a_nsp11_rx + ORF1a_nsp11.get_width()/2.0
    ORF1a_nsp11_cy = ORF1a_nsp11_ry + ORF1a_nsp11.get_height()/2.0

    # ORF1b
    ORF1b = patches.Rectangle((133.8 * 1.4, variant[1]), width=80.88 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1b)

    ORF1b_rx, ORF1b_ry = ORF1b.get_xy()
    ORF1b_cx = ORF1b_rx + ORF1b.get_width()/2.0
    ORF1b_cy = ORF1b_ry + ORF1b.get_height()/2.0

    ORF1a_nsp12 = patches.Rectangle((133.8 * 1.4, variant[1]), width=27.69 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1a_nsp12)

    ORF1a_nsp12_rx, ORF1a_nsp12_ry = ORF1a_nsp12.get_xy()
    ORF1a_nsp12_cx = ORF1a_nsp12_rx + ORF1a_nsp12.get_width()/2.0
    ORF1a_nsp12_cy = ORF1a_nsp12_ry + ORF1a_nsp12.get_height()/2.0

    ORF1a_nsp13 = patches.Rectangle((161.49 * 1.4, variant[1]), width=18.03 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1a_nsp13)

    ORF1a_nsp13_rx, ORF1a_nsp13_ry = ORF1a_nsp13.get_xy()
    ORF1a_nsp13_cx = ORF1a_nsp13_rx + ORF1a_nsp13.get_width()/2.0
    ORF1a_nsp13_cy = ORF1a_nsp13_ry + ORF1a_nsp13.get_height()/2.0

    ORF1a_nsp14 = patches.Rectangle((179.52 * 1.4, variant[1]), width=15.81 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1a_nsp14)

    ORF1a_nsp14_rx, ORF1a_nsp14_ry = ORF1a_nsp14.get_xy()
    ORF1a_nsp14_cx = ORF1a_nsp14_rx + ORF1a_nsp14.get_width()/2.0
    ORF1a_nsp14_cy = ORF1a_nsp14_ry + ORF1a_nsp14.get_height()/2.0

    ORF1a_nsp15 = patches.Rectangle((195.33 * 1.4, variant[1]), width=10.38 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1a_nsp15)

    ORF1a_nsp15_rx, ORF1a_nsp15_ry = ORF1a_nsp15.get_xy()
    ORF1a_nsp15_cx = ORF1a_nsp15_rx + ORF1a_nsp15.get_width()/2.0
    ORF1a_nsp15_cy = ORF1a_nsp15_ry + ORF1a_nsp15.get_height()/2.0

    ORF1a_nsp16 = patches.Rectangle((205.71 * 1.4, variant[1]), width=8.94 * 1.4, facecolor="#dce4f3", **kwargs)
    plt.gca().add_patch(ORF1a_nsp16)

    ORF1a_nsp16_rx, ORF1a_nsp16_ry = ORF1a_nsp16.get_xy()
    ORF1a_nsp16_cx = ORF1a_nsp16_rx + ORF1a_nsp16.get_width()/2.0
    ORF1a_nsp16_cy = ORF1a_nsp16_ry + ORF1a_nsp16.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((214.65 * 1.4, variant[1]), width=0.07 * 1.4, facecolor="#e0dcdb", **kwargs))

    # S
    S = patches.Rectangle((214.72 * 1.4, variant[1]), width=38.22 * 1.4, facecolor="#f7cfad", **kwargs)
    plt.gca().add_patch(S)

    S_rx, S_ry = S.get_xy()
    S_cx = S_rx + S.get_width()/2.0
    S_cy = S_ry + S.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((252.94 * 1.4, variant[1]), width=0.08 * 1.4, facecolor="#e0dcdb", **kwargs))

    # ORF3a
    ORF3a = patches.Rectangle((253.02 * 1.4, variant[1]), width=8.28 * 1.4, facecolor="#f1b9ba", **kwargs)
    plt.gca().add_patch(ORF3a)

    ORF3a_rx, ORF3a_ry = ORF3a.get_xy()
    ORF3a_cx = ORF3a_rx + ORF3a.get_width()/2.0
    ORF3a_cy = ORF3a_ry + ORF3a.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((261.3 * 1.4, variant[1]), width=0.24 * 1.4, facecolor="#e0dcdb", **kwargs))
    
    # E
    E = patches.Rectangle((261.54 * 1.4, variant[1]), width=2.28 * 1.4, facecolor="#d5edf2", **kwargs)
    plt.gca().add_patch(E)

    E_rx, E_ry = E.get_xy()
    E_cx = E_rx + E.get_width()/2.0
    E_cy = E_ry + E.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((263.82 * 1.4, variant[1]), width=0.5 * 1.4, facecolor="#e0dcdb", **kwargs))

    # M
    M = patches.Rectangle((264.32 * 1.4, variant[1]), width=6.69 * 1.4, facecolor="#bbd6b6", **kwargs)
    plt.gca().add_patch(M)

    M_rx, M_ry = M.get_xy()
    M_cx = M_rx + M.get_width()/2.0
    M_cy = M_ry + M.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((271.01 * 1.4, variant[1]), width=0.1 * 1.4, facecolor="#e0dcdb", **kwargs))

    # ORF6
    ORF6 = patches.Rectangle((271.11 * 1.4, variant[1]), width=1.86 * 1.4, facecolor="#f4e5b7", **kwargs)
    plt.gca().add_patch(ORF6)

    ORF6_rx, ORF6_ry = ORF6.get_xy()
    ORF6_cx = ORF6_rx + ORF6.get_width()/2.0
    ORF6_cy = ORF6_ry + ORF6.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((272.97 * 1.4, variant[1]), width=0.06 * 1.4, facecolor="#e0dcdb", **kwargs))

    # ORF7a
    ORF7a = patches.Rectangle((273.03 * 1.4, variant[1]), width=3.66 * 1.4, facecolor="#d2bee3", **kwargs)
    plt.gca().add_patch(ORF7a)

    ORF7a_rx, ORF7a_ry = ORF7a.get_xy()
    ORF7a_cx = ORF7a_rx + ORF7a.get_width()/2.0
    ORF7a_cy = ORF7a_ry + ORF7a.get_height()/2.0

    # ORF7b
    ORF7b = patches.Rectangle((276.69 * 1.4, variant[1]), width=1.32 * 1.4, facecolor="#fcd5d9", **kwargs)
    plt.gca().add_patch(ORF7b)

    ORF7b_rx, ORF7b_ry = ORF7b.get_xy()
    ORF7b_cx = ORF7b_rx + ORF7b.get_width()/2.0
    ORF7b_cy = ORF7b_ry + ORF7b.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((278.01 * 1.4, variant[1]), width=0.06 * 1.4, facecolor="#e0dcdb", **kwargs))
    
    # ORF8
    ORF8 = patches.Rectangle((278.07 * 1.4, variant[1]), width=3.66 * 1.4, facecolor="#cdb7b4", **kwargs)
    plt.gca().add_patch(ORF8)

    ORF8_rx, ORF8_ry = ORF8.get_xy()
    ORF8_cx = ORF8_rx + ORF8.get_width()/2.0
    ORF8_cy = ORF8_ry + ORF8.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((281.73 * 1.4, variant[1]), width=0.14 * 1.4, facecolor="#e0dcdb", **kwargs))

    # N
    N = patches.Rectangle((281.87 * 1.4, variant[1]), width=12.6 * 1.4, facecolor="#b7b7b8", **kwargs)
    plt.gca().add_patch(N)

    N_rx, N_ry = N.get_xy()
    N_cx = N_rx + N.get_width()/2.0
    N_cy = N_ry + N.get_height()/2.0

    # ORF9b
    ORF9b = patches.Rectangle((281.97 * 1.4, variant[1]), width=2.94 * 1.4, facecolor="#9fc9f3", **kwargs)
    plt.gca().add_patch(ORF9b)

    ORF9b_rx, ORF9b_ry = ORF9b.get_xy()
    ORF9b_cx = ORF9b_rx + ORF9b.get_width()/2.0
    ORF9b_cy = ORF9b_ry + ORF9b.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((294.47 * 1.4, variant[1]), width=0.23 * 1.4, facecolor="#e0dcdb", **kwargs))
    
    # ORF10
    ORF10 = patches.Rectangle((294.7 * 1.4, variant[1]), width=1.17 * 1.4, facecolor="#dfe1ad", **kwargs)
    plt.gca().add_patch(ORF10)

    ORF10_rx, ORF10_ry = ORF10.get_xy()
    ORF10_cx = ORF10_rx + ORF10.get_width()/2.0
    ORF10_cy = ORF10_ry + ORF10.get_height()/2.0

    # n/a
    plt.gca().add_patch(patches.Rectangle((295.87 * 1.4, variant[1]), width=2.29 * 1.4, facecolor="#e0dcdb", **kwargs))

    # Nonsynonymous
    ORF1a = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF1a")]

    for index, row in ORF1a.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF1a")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF1a")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#a7c5df"), color="#a7c5df", fontproperties=font_properties)

    ORF1b = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF1b")]

    for index, row in ORF1b.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF1b")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF1b")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#dce4f3"), color="#dce4f3", fontproperties=font_properties)

    S = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "S")]

    for index, row in S.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("S")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("S")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#f7cfad"), color="#f7cfad", fontproperties=font_properties)

    ORF3a = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF3a")]

    for index, row in ORF3a.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF3a")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF3a")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#f1b9ba"), color="#f1b9ba", fontproperties=font_properties)

    E = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "E")]

    for index, row in E.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("E")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("E")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#d5edf2"), color="#d5edf2", fontproperties=font_properties)

    M = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "M")]

    for index, row in M.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("M")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("M")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#bbd6b6"), color="#bbd6b6", fontproperties=font_properties)

    # There's no ORF6

    ORF7a = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF7a")]

    for index, row in ORF7a.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF7a")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF7a")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#d2bee3"), color="#d2bee3", fontproperties=font_properties)

    ORF7b = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF7b")]

    for index, row in ORF7b.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF7b")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF7b")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#fcd5d9"), color="#fcd5d9", fontproperties=font_properties)

    ORF8 = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF8")]

    for index, row in ORF8.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("ORF8")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("ORF8")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#cdb7b4"), color="#cdb7b4", fontproperties=font_properties)

    N = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "N")]

    for index, row in N.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("N")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("N")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#b7b7b8"), color="#b7b7b8", fontproperties=font_properties)

    ORF9b = df.loc[df["gene"].notnull()].loc[np.logical_and(df["classification"] == variant[0], df["gene"] == "ORF9b")]

    for index, row in ORF9b.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("N")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("N")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="#9fc9f3"), color="#9fc9f3", fontproperties=font_properties)

    # Synonymous
    synonymous = df.loc[df["gene"].isnull()].loc[df["classification"] == variant[0]]
    for index, row in synonymous.iterrows():
        amino_acid, position = row["amino_acid"], row["position"]

        plt.annotate(amino_acid, rotation=90, xy=(formulas.get("start")(position) * 1.4, variant[1] + .5), xytext=(formulas.get("start")(position) * 1.4, variant[1] + 2), xycoords="data", 
                    fontsize=8, ha="center", va="top",
                    arrowprops=dict(arrowstyle="-", lw=2.0, color="black"), color="black", fontproperties=font_properties)

width = 8
custom = [Line2D([0], [0], color="#a7c5df", lw=width),
         Line2D([0], [0], color="#dce4f3", lw=width),
         Line2D([0], [0], color="#f7cfad", lw=width),
         Line2D([0], [0], color="#f1b9ba", lw=width),
         Line2D([0], [0], color="#d5edf2", lw=width),
         Line2D([0], [0], color="#bbd6b6", lw=width),
         Line2D([0], [0], color="#f4e5b7", lw=width),
         Line2D([0], [0], color="#d2bee3", lw=width),
         Line2D([0], [0], color="#fcd5d9", lw=width),
         Line2D([0], [0], color="#cdb7b4", lw=width),
         Line2D([0], [0], color="#b7b7b8", lw=width),
         Line2D([0], [0], color="#9fc9f3", lw=width),
         Line2D([0], [0], color="#dfe1ad", lw=width)]

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

ax.text(-15.8 * 1.4, 16.1, "Alpha", fontproperties=font_properties, fontsize=16)
ax.text(-13.7 * 1.4, 12.1, "Beta", fontproperties=font_properties, fontsize=16)
ax.text(-14.8 * 1.4, 8.1, "Delta", fontproperties=font_properties, fontsize=16)
ax.text(-20.2 * 1.4, 4.1, "Gamma", fontproperties=font_properties, fontsize=16)
ax.text(-22 * 1.4, .1, "Omicron", fontproperties=font_properties, fontsize=16)

ax.legend(custom, ["ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF10"], title="Genes", borderpad=2, bbox_to_anchor=(-.01, 1.0), loc="upper right")
plt.tight_layout()
# plt.savefig("plot_2.pdf", transparent=True)
# subprocess.Popen(["open", "plot.pdf"])
plt.show()

