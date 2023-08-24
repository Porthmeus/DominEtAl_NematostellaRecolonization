# Porthmeus
# 09.06.20

rule plotVenn_custom:
    input:
        venn = "DESeq2/VennMatrix.csv",
        vennAnno = "DESeq2/VennAnnoTable.csv"
    output:
        venn = "Plots/VennCustom.svg",
        geneList = "Data/vennAreas_GF_vs_bA_bJ_bL.txt",
        geneTable = "Data/DEG_vennAnno.csv"
    conda: "envs/R.yaml"
    log: "logs/plotVenn_custom.log"
    script: "../scripts/plotVenn_custom.R"
    
