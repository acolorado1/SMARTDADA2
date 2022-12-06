rule all: 
    input: 
        "output/parameter_info_LOZ_nano_fixed.tsv",
        "output/sumEEInfo_LOZ_nano.tsv",
        "output/ScatterReadLengthByAvgEE.png",
        "output/HeatmapIndexValueByAvgEE.png",
        "output/HistogramRetainedReadCount.png"

rule create_tsvs:
    input: 
        "LOZ-CSU-Nano_S1_L001_R1_001.fastq"
    output: 
        TrimInfo = "output/parameter_info_LOZ_nano_fixed.tsv",
        SumEEInfo = "output/sumEEInfo_LOZ_nano.tsv"
    params: 
        threshold = 30, 
        o_mtp = 0.1, 
        a_mtp = 0.2
    shell: 
        "python " + "smartdada2/main.py --fq {input} --t_of {output.TrimInfo} --th {params.threshold} --o_mtp {params.o_mtp} --a_mtp {params.a_mtp} --EE_of {output.SumEEInfo}"

rule get_plots: 
    input: 
        TrimInfo = "output/parameter_info_LOZ_nano_fixed.tsv", 
        sumEEInfo = "output/sumEEInfo_LOZ_nano.tsv"
    output:
        "output/ScatterReadLengthByAvgEE.png",
        "output/HeatmapIndexValueByAvgEE.png",
        "output/HistogramRetainedReadCount.png"
    shell:
        "Rscript " + "R_scripts/vizualize.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo}"
