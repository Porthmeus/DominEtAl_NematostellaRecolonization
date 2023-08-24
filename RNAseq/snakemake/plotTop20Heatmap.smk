# Porthmeus
# 24.10.19
rule plotTop20Heatmap:
    input:
        res = "DESeq2/{comp}.csv",
        deseq = "DESeq2/DeseqMat.RData"
    output:
        top = report("Plots/resTop40Heat_{comp}.svg",
            caption = "../report/plotResTop40Heat.rst",
            category = "Results"),
        up = report("Plots/resTop20UpHeat_{comp}.svg",
            caption = "../report/plotResTop20UpHeat.rst",
            category = "Results"),
        down = report("Plots/resTop20DownHeat_{comp}.svg",
            caption = "../report/plotResTop20DownHeat.rst",
            category = "Results"),
    conda:
        "envs/R.yaml"
    log:
        "logs/resultsTop20Heat_{comp}.log"
    script:
        "../scripts/plotTop20Heatmap.R"
