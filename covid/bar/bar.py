from matplotlib import pyplot as plt
import matplotlib
import json

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

data = {
   "alpha" : {
        # ORF1a
        "nsp1" : 0,
        "nsp2" : 0,
        "nsp3" : 0,
        "nsp4" : 0,
        "nsp5" : 0,
        "nsp6" : 0,
        "nsp7" : 0,
        "nsp8" : 0,
        "nsp9" : 0,
        "nsp10" : 0,
        "nsp11" : 0,
        # ORF1b
        "nsp12" : 0,
        "nsp13" : 0,
        "nsp14" : 0,
        "nsp15" : 0,
        "nsp16" : 0,
        "S" : 0, 
        "ORF3a" : 0,
        "E" : 0,
        "M" : 0,
        "ORF6" : 0,
        "ORF7a" : 0,
        "ORF7b" : 0,
        "ORF8" : 0,
        "N" : 0,
        "ORF9b" : 0, 
        "ORF10" : 0
   },
   "beta" : {
        # ORF1a
        "nsp1" : 0,
        "nsp2" : 0,
        "nsp3" : 0,
        "nsp4" : 0,
        "nsp5" : 0,
        "nsp6" : 0,
        "nsp7" : 0,
        "nsp8" : 0,
        "nsp9" : 0,
        "nsp10" : 0,
        "nsp11" : 0,
        # ORF1b
        "nsp12" : 0,
        "nsp13" : 0,
        "nsp14" : 0,
        "nsp15" : 0,
        "nsp16" : 0,
        "S" : 0, 
        "ORF3a" : 0,
        "E" : 0,
        "M" : 0,
        "ORF6" : 0,
        "ORF7a" : 0,
        "ORF7b" : 0,
        "ORF8" : 0,
        "N" : 0,
        "ORF9b" : 0, 
        "ORF10" : 0
   },
   "delta" : {
        # ORF1a
        "nsp1" : 0,
        "nsp2" : 0,
        "nsp3" : 0,
        "nsp4" : 0,
        "nsp5" : 0,
        "nsp6" : 0,
        "nsp7" : 0,
        "nsp8" : 0,
        "nsp9" : 0,
        "nsp10" : 0,
        "nsp11" : 0,
        # ORF1b
        "nsp12" : 0,
        "nsp13" : 0,
        "nsp14" : 0,
        "nsp15" : 0,
        "nsp16" : 0,
        "S" : 0, 
        "ORF3a" : 0,
        "E" : 0,
        "M" : 0,
        "ORF6" : 0,
        "ORF7a" : 0,
        "ORF7b" : 0,
        "ORF8" : 0,
        "N" : 0,
        "ORF9b" : 0, 
        "ORF10" : 0
   },
   "gamma" : {
        # ORF1a
        "nsp1" : 0,
        "nsp2" : 0,
        "nsp3" : 0,
        "nsp4" : 0,
        "nsp5" : 0,
        "nsp6" : 0,
        "nsp7" : 0,
        "nsp8" : 0,
        "nsp9" : 0,
        "nsp10" : 0,
        "nsp11" : 0,
        # ORF1b
        "nsp12" : 0,
        "nsp13" : 0,
        "nsp14" : 0,
        "nsp15" : 0,
        "nsp16" : 0,
        "S" : 0, 
        "ORF3a" : 0,
        "E" : 0,
        "M" : 0,
        "ORF6" : 0,
        "ORF7a" : 0,
        "ORF7b" : 0,
        "ORF8" : 0,
        "N" : 0,
        "ORF9b" : 0, 
        "ORF10" : 0
   },
   "omicron" : {
        # ORF1a
        "nsp1" : 0,
        "nsp2" : 0,
        "nsp3" : 0,
        "nsp4" : 0,
        "nsp5" : 0,
        "nsp6" : 0,
        "nsp7" : 0,
        "nsp8" : 0,
        "nsp9" : 0,
        "nsp10" : 0,
        "nsp11" : 0,
        # ORF1b
        "nsp12" : 0,
        "nsp13" : 0,
        "nsp14" : 0,
        "nsp15" : 0,
        "nsp16" : 0,
        "S" : 0, 
        "ORF3a" : 0,
        "E" : 0,
        "M" : 0,
        "ORF6" : 0,
        "ORF7a" : 0,
        "ORF7b" : 0,
        "ORF8" : 0,
        "N" : 0,
        "ORF9b" : 0, 
        "ORF10" : 0
   }
}

variants = ["Alpha", "Beta", "Delta", "Gamma", "Omicron"]
synonymous = [7, 4, 7, 10, 10]
nonsynonymous = [24, 20, 29, 24, 62]

plt.barh(variants, synonymous, color="#682C8B")
plt.barh(variants, nonsynonymous, left=synonymous, color="#BEBDBE")

plt.show()