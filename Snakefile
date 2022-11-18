rule all: 
    input: 
        "Report.pdf"

rule create_tsv:
    input: 
        "data/SRR1591840.fastq"
    output: 
        "data/parameter_info.tsv"
    params: 
        threshold = 30, 
        o_mtp = 0.1, 
        a_mtp = 0.2, 
        maxEE = 2.0
    shell: 
        "python " + "smartdada2/main.py --fq {input} --of {output} --th {params.threshold} --o_mtp {params.o_mtp} --a_mtp {params.a_mtp} --mEE {params.maxEE}"

rule get_plots: 
    input: 
        "data/parameter_info.tsv"
    output:
        "plots/ReadLengthByAvgEE.png",
        "plots/ReadCountOverMaxEE.png"
    shell:
        "Rscript " + "R_scripts/viz_SizeByAvgEE.R --tsv_file {input}"

rule update_report: 
    input: 
        "plots/ReadLengthByAvgEE.png",
        "plots/ReadCountOverMaxEE.png"
    output: 
        "Report.pdf"
    shell:
        "pandoc " + "Report.md -o {output}"
