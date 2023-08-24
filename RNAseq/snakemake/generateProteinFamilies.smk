# Porthmeus
# 18.03.20

rule generateProteinFamilies:
    input:
        ref = "proteinFamilies/Blast/nematostellaTx_vienna.pep",
        evalu = "proteinFamilies/Evaluation/Evaluation.csv",
        rmd = "proteinFamilies/NatrixDocumentation_ProteinFamilies.Rmd"
    output:
        dump = expand("proteinFamilies/TribeMCL/dump.seq.mci.I{inflation}", inflation = ["14","20","40","60","100"]),
        diagHist = "proteinFamilies/Plots/diagHist.svg",
        diagBar = "proteinFamilies/Plots/diagBar.svg",
        html = report("proteinFamilies/NatrixDocumentation_ProteinFamilies.html",
                caption = "../report/NatrixDocumentation_ProteinFamilies.rst",
                category = "ProteinFamily")
    log:
        log="logs/NatrixDocumentation_ProteinFamilies.log"
    conda:
        "envs/proteinFamilies.yaml"
    shell:
        '''R -e "rmarkdown::render('{input.rmd}')" 2&> {log}'''
