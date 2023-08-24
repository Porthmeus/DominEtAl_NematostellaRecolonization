# 28.02.18
# Porthmeus
import pandas as pd


df = {"PF" : [],"Gene" : []}
i=0
with open(snakemake.input[0], "r") as infile:
    for line in infile:
        i=i+1
        PF = "PF"+str(i)
        l=line.strip("\n").split("\t")
        for gene in l:
            g = gene.split("::")[1]
            df["PF"].append(PF)
            df["Gene"].append(g)

df = pd.DataFrame(df, columns = df.keys()).drop_duplicates(subset="Gene")
df.to_csv(snakemake.output[0], index=False)
