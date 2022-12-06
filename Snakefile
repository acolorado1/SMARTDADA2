rule all: 
    input: 
        "output/parameter_info_LOZ_nano_fixed.tsv",
        "output/sumEEInfo_LOZ_nano.tsv"

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

'''rule get_plots: 
    input: 
        "data/parameter_info_LOZ_nano_fixed.tsv"
    output:
        "output/ReadLengthByAvgEE.png",
        "output/ReadCountOverMaxEE.png",
        "output/t_len_by_AvgEE.png",
        "output/t_len_by_ReadsOverMaxEE.png"
    shell:
        "Rscript " + "R_scripts/vizualize.R --tsv_file {input}"'''
