# Porthmeus
# 24.10.19
rule plotDiagCorHeatmap:
    input:
        deseq = "DESeq2/DeseqMat.RData"
    output:
        wb0 =   report(  "Plots/DiagCorHeatmap_wb0.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wb1 =   report(  "Plots/DiagCorHeatmap_wb1.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wb3 =   report(  "Plots/DiagCorHeatmap_wb3.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wb10 =  report( "Plots/DiagCorHeatmap_wb10.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wob0 =  report( "Plots/DiagCorHeatmap_wob0.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wob1 =  report( "Plots/DiagCorHeatmap_wob1.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wob3 =  report( "Plots/DiagCorHeatmap_wob3.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
        wob10 = report("Plots/DiagCorHeatmap_wob10.svg",
            caption = "../report/plotDiagCorHeatmap.rst",
            category = "Diagnostics"),
    log:
        "logs/DiagCorHeatmapsPlots.log"
    conda:
        "envs/plotDiagCorHeatmaps.yaml"
    script:
        "../scripts/plotDiagCorHeatmap.R"
