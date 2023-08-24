# Porthmeus
# 23.10.19

rule tximport:
    input:
        meta = "metaAll.csv"
    output:
        "DESeq2/DeseqMat.RData"
    params:
        files = "mapping/"
    log:
        "logs/tximport.log"
    conda:
        "envs/R.yaml"
    script:
        "../scripts/tximport.R"
