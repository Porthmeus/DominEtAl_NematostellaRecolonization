rule testSubsetKEGGEnrichment:
    input:
        deseq = "DESeq2/DeseqMat.RData",
        results = "Data/subset_{compSub}.txt",
        cluster = "Annotation/Cluster2KeggLongWhite.csv"
    output:
        hyGeoPlot = report("KEGG/Enrichment/SubsetPlot_{compSub}_KEGG_hyGeo.svg",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
        hyGeoCSV = report("KEGG/Enrichment/SubsetData_{compSub}_KEGG_hyGeo.csv",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
        hyGeoRData = report("KEGG/Enrichment/SubsetData_{compSub}_KEGG_hyGeo.RData",
            caption = "../report/hyGeoKEGG.rst",
            category = "KEGG"),
    conda:
        "envs/R.yaml"
    log:
        "logs/{compSub}_test_KEGG_EnrichmentSubset.log"
    script:
        "../scripts/testEnrichmentSubset.R"
        

