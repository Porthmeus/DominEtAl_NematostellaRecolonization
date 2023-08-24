# Porthmeus
# 11.06.20

# split the gene list of the custom kegg plots
import re

rule subsetGenelist:
    input: 
        lst = "Data/vennAreas_GF_vs_bA_bJ_bL.txt"
    output: 
        expand("Data/subset_{comp}.txt", 
                comp = [line.strip().replace("# ","") for line in open("Data/vennAreas_GF_vs_bA_bJ_bL.txt","r") if line.startswith("#")])
    log: "logs/subsetGeneList.log"
    run:
        def grep(pattern, lst):
            pt = re.compile(".*(" + pattern + ").*")
            return [m.group(0) for l in lst for m in [pt.search(l)] if m]
        fl = None
        with open(input[0], "r") as inp:
            for line in inp:
                if line.startswith("#"):
                    fl = grep(line.strip().replace("# ",""), output)[0]

                else:
                    if fl != None:
                        with open(fl, "a") as out:
                            out.write(line)
