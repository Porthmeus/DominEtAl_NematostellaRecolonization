# Porthmeus
# 06.02.20

rule ResultsToKegg:
    input:
        results = expand("DESeq2/{comp}.csv", comp = constructComparison(comparisons))
    output:
        kegg = report("DESeq2/KEGGMapping.txt",
            caption = "../report/KEGGMapping.rst",
            category = "Results")
    conda:
        "envs/R.yaml"
    log:
        "logs/ResultsToKegg.log"
    script:
        "../scripts/ResultsToKegg.R"
