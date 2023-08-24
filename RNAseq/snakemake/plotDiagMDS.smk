# Porthmeus
# 24.10.19
rule plotDiagMDS:
    input:
        deseq = "DESeq2/DeseqMat.RData"
    output:
        mds = report("Plots/DiagMDS.svg",
            caption = "../report/plotDiagMDS.rst",
            category = "Diagnostics"),
        stress = report("Plots/DiagMDSStress.svg",
            caption = "../report/plotDiagMDSStress.rst",
            category = "Diagnostics"),
        csv = report("Data/dataDiagMDS.csv",
            caption = "../report/dataDiagMDS.rst",
            category = "Diagnostics")
    log:
        "logs/DiagMDSPlots.log"
    conda:
        "envs/R.yaml"
    script:
        "../scripts/plotDiagMDS2.R"
