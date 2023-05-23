FASTQ_FILE="smartdada2/testing/test_data/LOZ_Nano_Trunc.fastq"
PAIRED=True
SUBSAMPLE=500000

if PAIRED:
    inputs=["forward", "reverse"]
else:
    inputs=["forward"]

rule all:
    input:
        "forward/output/SMARTDADA2_InteractiveOutput.html",
        "reverse/output/SMARTDADA2_InteractiveOutput.html"

if PAIRED: 
    rule get_paried_files:
        input: 
            fq=FASTQ_FILE
        output:
            forward=inputs[0] + "/forward.fastq",
            rev=inputs[1] + "/reverse.fastq"
        shell:
            "python " + "smartdada2/GetPairedendFiles.py --fq {input.fq} --of_f {output.forward} --of_r {output.rev}"

rule create_TSVs:
    input:
        "{inputs}/{inputs}.fastq"
    output:
        TrimInfo = "{inputs}/output/TrimInfo.tsv",
        SumEEInfo = "{inputs}/output/sumEEInfo.tsv"
    params:
        threshold = 30,
        o_mtp = 0.1,
        a_mtp = 0.2
    shell:
        "python " + "smartdada2/main.py --fq {input} --t_of {output.TrimInfo} --th {params.threshold} --o_mtp {params.o_mtp} --a_mtp {params.a_mtp} --EE_of {output.SumEEInfo}"

rule get_plots:
    input:
        TrimInfo = "{inputs}/output/TrimInfo.tsv",
        sumEEInfo = "{inputs}/output/sumEEInfo.tsv",
        output_file = "{inputs}"
    output:
        "{inputs}/output/plots/ScatterReadLengthByAvgEE.png",
        "{inputs}/output/plots/HeatmapIndexValueByAvgEE.png",
        "{inputs}/output/plots/HistogramRetainedReadCount.png"
    shell:
        "Rscript " + "smartdada2/R_scripts/vizualize.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo} --output_file {input.output_file}"

rule get_interactive_plots:
    input:
        TrimInfo = "{inputs}/output/TrimInfo.tsv",
        sumEEInfo = "{inputs}/output/sumEEInfo.tsv",
        output_file = "{inputs}"
    output:
        "{inputs}/output/SMARTDADA2_InteractiveOutput.html"
    shell:
        "Rscript " + "smartdada2/R_scripts/wrapper_RunRender.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo} --output_file {input.output_file}"



