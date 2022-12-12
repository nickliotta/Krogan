import pandas as pd
import json
new = []

with open("data.json") as file:
    file = json.load(file)
    for classification in file:
        synonymous = file[classification].get("synonymous")
        for i in synonymous:
            new.append({
                "classification" : classification,
                "mutation" : "synonymous",
                "amino_acid" : i,
                "position" : int(''.join(filter(str.isdigit, i)))
            })

with open("data.json") as file:
    file = json.load(file)
    for classification in file:
        for gene in file[classification].get("nonsynonymous"):
            for amino in file[classification].get("nonsynonymous")[gene]:
                new.append({
                    "classification" : classification,
                    "gene" : gene,      
                    "mutation" : "nonsynonymous",
                    "amino_acid" : amino,
                    "position" : int(''.join(filter(str.isdigit, amino))),
                })

df = pd.DataFrame.from_dict(new) 
df.to_csv("data.csv", index=False, header=True)