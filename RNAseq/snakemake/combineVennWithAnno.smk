# Porthmeus
# 10.06.20

# combine the generated Venn Matrix with the Annotation

rule combineVennWithAnno:
    input:
        venn = "DESeq2/VennDirectionalMatrix.csv",
        anno = "Annotation/Contig2AnnoAggregated.csv"
    output:
        vennAnno = "DESeq2/VennAnnoTable.csv"
    log:"logs/combineVennWithAnno.log"
    conda: "../envs/R.yaml"
    script: "../scripts/combineVennWithAnno.R"
