# Porthmeus
# 24.10.19
rule plotDiagPCA:
    input:
        deseq = "DESeq2/DeseqMat.RData"
    output:
        PCA = report( "Plots/DiagPCA.svg", 
            caption = "../report/plotDiagPCA.rst",
            category = "Diagnostics"),
        csv = report("Data/dataDiagPCA.csv",
            caption = "../report/dataDiagPCA.rst",
            category = "Diagnostics")
    log:
        "logs/DiagPCAPlots.log"
    conda:
        "envs/R.yaml"
    script:
        "../scripts/plotDiagPCA.R"
