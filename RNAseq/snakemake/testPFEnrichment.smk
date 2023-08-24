# Porthmeus
# 20.03.20

rule testPFEnrichment:
    input:
        results = "DESeq2/{comp}.RData",
        cluster = "proteinFamilies/Enrichment/ProteinFamilies.csv"
    output:
        hyGeoPlot = report("proteinFamilies/Enrichment/plot{comp}_PF_hyGeo.svg",
            caption = "../report/hyGeoPF.rst",
            category = "ProteinFamily"),
        hyGeoCSV = report("proteinFamilies/Enrichment/data{comp}_PF_hyGeo.csv",
            caption = "../report/hyGeoPF.rst",
            category = "ProteinFamily"),
        hyGeoRData = report("proteinFamilies/Enrichment/data{comp}_PF_hyGeo.RData",
            caption = "../report/hyGeoPF.rst",
            category = "ProteinFamily"),
        gseaPlot = report("proteinFamilies/Enrichment/plot{comp}_PF_gsea.svg",
            caption = "../report/gseaPF.rst",
            category = "ProteinFamily"),
        gseaCSV = report("proteinFamilies/Enrichment/data{comp}_PF_gsea.csv",
            caption = "../report/gseaPF.rst",
            category = "ProteinFamily"),
        gseaRData = report("proteinFamilies/Enrichment/data{comp}_PF_gsea.RData",
            caption = "../report/gseaPF.rst",
            category = "ProteinFamily")
    conda:
        "envs/R.yaml"
    log:
        "logs/{comp}_test_KEGG_Enrichment.log"
    script:
        "../scripts/testEnrichment.R"
        
