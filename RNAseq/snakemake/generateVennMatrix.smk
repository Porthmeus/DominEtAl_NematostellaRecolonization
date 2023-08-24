# Porthmeus
# 09.06.20

# generate a venn matrix for all comparisons
rule generateVennMatrix:
    input:
        expand("DESeq2/{comp}.RData", comp = constructComparison(comparisons))
    output:
        venn = "DESeq2/VennMatrix.csv",
        vennDir = "DESeq2/VennDirectionalMatrix.csv"
    log: "logs/generateVennMatrix.log"
    conda: "envs/R.yaml"
    script: "../scripts/generateVennMatrix.R"
