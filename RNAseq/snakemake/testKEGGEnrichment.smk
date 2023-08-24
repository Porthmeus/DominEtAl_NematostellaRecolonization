# Porthmeus
# 20.03.20

rule testKEGGEnrichment:
    input:
        results = "DESeq2/{comp}.RData",
        cluster = "Annotation/Cluster2KeggLongWhite.csv"
    output:
        hyGeoPlot = report("KEGG/Enrichment/plot{comp}_KEGG_hyGeo.svg",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
        hyGeoCSV = report("KEGG/Enrichment/data{comp}_KEGG_hyGeo.csv",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
        hyGeoRData = report("KEGG/Enrichment/data{comp}_KEGG_hyGeo.RData",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
        gseaPlot = report("KEGG/Enrichment/plot{comp}_KEGG_gsea.svg",
            caption = "../report/gseaKEGG.rst",
            category = "KEGG"),
        gseaCSV = report("KEGG/Enrichment/data{comp}_KEGG_gsea.csv",
            caption = "../report/gseaKEGG.rst",
            category = "KEGG"),
        gseaRData = report("KEGG/Enrichment/data{comp}_KEGG_gsea.RData",
            caption = "../report/gseaKEGG.rst",
            category = "KEGG")
    conda:
        "envs/R.yaml"
    log:
        "logs/{comp}_test_KEGG_Enrichment.log"
    script:
        "../scripts/testEnrichment.R"
        
