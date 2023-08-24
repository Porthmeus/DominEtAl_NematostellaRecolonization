# Porthmeus
# 24.10.19
rule DESeq2_results:
    input:
        deseq = "DESeq2/DeseqMat.RData",
        anno = "Annotation/Cluster2Kegg.csv"
    output:
        results = report("DESeq2/{comp}.csv",
            caption = "../report/resDESeq2results.rst",
            category = "Results"),
        resultsR = report("DESeq2/{comp}.RData",
            caption = "../report/resDESeq2results.rst",
            category = "Results")
    conda:
        "envs/R.yaml"
    log:
        "logs/DESeq2_results_{comp}.log"
    script:
        "../scripts/DESeq2_results2.R"
