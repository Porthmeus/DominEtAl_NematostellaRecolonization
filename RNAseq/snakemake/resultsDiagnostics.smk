# Porthmeus
# 24.10.19
rule resultsDiagnostics:
    input:
        "DESeq2/{comp}.RData",
    output:
        MA = report("Plots/resDiagMAPlot_{comp}.jpg",
            caption = "../report/plotResMAPlot.rst",
            category = "Results"),
        Hists = report("Plots/resDiagHistPvalsFCbaseMeanPlot_{comp}.svg",
            caption = "../report/plotResDiagHistPvalsFCbaseMeanPlot.rst",
            category = "Results"),
        SmallP = report("Plots/resDiagSmallPToCountsPlot_{comp}.svg",
            caption = "../report/plotDiagAmallPToCountsPlot.rst",
            category = "Results")
    conda:
        "envs/R.yaml"
    log:
        "logs/resultsDiagnostics_{comp}.log"
    script:
        "../scripts/resultsDiagnostics.R"
