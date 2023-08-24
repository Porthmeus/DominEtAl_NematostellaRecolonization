# Porthmeus
# 23.10.19

#report:"report/workflow.rst"

# define here the comparisons which should be done throughout the analysis. This
# can be quite a few as experimental design increases. The key of the dictionary
# is the name of the column in the meta data sheet, while the values are the
# actual comparisons. The keys must be part of the design table in the tximport
# script. It is not allowed to use "_" as delimiter here, otherwise an ERROR
# will occur.
comparisons = {"Condition":[[["GF"],["wt"]],
                [["GF"],["bA.A"]],
                [["GF"],["bJ.A"]],
                [["GF"],["bE.A"]],
                [["bA.A"],["wt"]],
                [["bJ.A"],["wt"]],
                [["bE.A"],["wt"]],
                [["bE.A","bJ.A","bA.A"],["wt"]],
                [["bE.A","bJ.A","bA.A"],["GF"]],
                [["bE.A","bJ.A","bA.A"],["GF","wt"]]
                ]}
# function to construct the comparisons which can then be parsed to the
# different functions throughout the analysis
def constructComparison(comparison):
    comp = []
    for key in comparison.keys():
        for vals in comparison[key]:
            if len(vals[0]) > 1:
                vals[0] = ["_".join(vals[0])]
            if len(vals[1]) > 1:
                vals[1] = ["_".join(vals[1])]
            txt = key + "_" + vals[0][0] + "_vs_" + vals[1][0]
            comp.append(txt)
    return (comp)

rule all:
    input:
        expand([
        "Plots/DiagPCA.svg",
        "Plots/DiagMDS.svg",
        "Plots/DiagCorHeatmap_{wb}{thrld}.svg",
        "Plots/resTop40Heat_{comp}.svg",
        "Plots/resTop20UpHeat_{comp}.svg",
        "Plots/resTop20DownHeat_{comp}.svg",
        "Plots/resDiagMAPlot_{comp}.jpg",
        "Plots/resDiagHistPvalsFCbaseMeanPlot_{comp}.svg",
        "Plots/resDiagSmallPToCountsPlot_{comp}.svg",
        "DESeq2/KEGGMapping.txt",
        "KEGG/Enrichment/plot{comp}_KEGG_gsea.svg",
        "KEGG/Enrichment/data{comp}_KEGG_gsea.csv",
        "KEGG/Enrichment/data{comp}_KEGG_gsea.RData",
        "KEGG/Enrichment/plot{comp}_KEGG_hyGeo.svg",
        "KEGG/Enrichment/data{comp}_KEGG_hyGeo.csv",
        "KEGG/Enrichment/data{comp}_KEGG_hyGeo.RData",
        "DESeq2/VennMatrix.csv",
        "Plots/VennCustom.svg",
        "DESeq2/VennAnnoTable.csv",
        "KEGG/Enrichment/SubsetPlot_{compSub}_KEGG_hyGeo.svg",
        ],  wb = ["wb","wob"],
            thrld = [0,1,3,10],
            comp = constructComparison(comparisons), 
            compSub = [
                "overlap",
                "Condition_GF_vs_bA.A",
                "Condition_GF_vs_bJ.A",
                "Condition_GF_vs_bE.A"]
            )


include: "snakemake/tximport.smk"
include: "snakemake/plotDiagPCA.smk"
include: "snakemake/plotDiagMDS.smk"
include: "snakemake/plotDiagCorHeatmap.smk"
include: "snakemake/DESeq2_results2.smk"
include: "snakemake/resultsDiagnostics.smk"
include: "snakemake/plotTop20Heatmap.smk"
include: "snakemake/ResultsToKegg.smk"
include: "snakemake/generateProteinFamilies.smk"
include: "snakemake/reformatProteinFamilies.smk"
include: "snakemake/testPFEnrichment.smk"
include: "snakemake/testKEGGEnrichment.smk"
include: "snakemake/generateVennMatrix.smk"
include: "snakemake/plotVenn_custom.smk"
include: "snakemake/combineVennWithAnno.smk"
include: "snakemake/subsetGenelists.smk"
include: "snakemake/testSubsetKegg.smk"
