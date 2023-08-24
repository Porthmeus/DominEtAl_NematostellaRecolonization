# Porthmeus
# 17.03.20

rule reformatProteinFamilies:
    input:
        "proteinFamilies/TribeMCL/dump.seq.mci.I40"
    output:
        report("proteinFamilies/Enrichment/ProteinFamilies.csv",
            caption = "../report/reformatProteinFamilies.rst",
            category = "ProteinFamily")
    log:
        log = "../logs/reformatProteinFamilies.log"
    conda:
        "envs/python.yaml"
    script:
        "../scripts/reformatProteinFamilies.py"
