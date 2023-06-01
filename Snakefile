FASTQ_FILE = "../DADA2ParameterExploration/smartdada2/testing/test_data/LOZ_Nano_Trunc.fastq"
PAIRED = True
SUBSAMPLE = 2
AVG_Q_SCORE = 30
OBV_TRIMMING_MAX = 0.1
MAX_TRIMMING = 0.2

if PAIRED:
    INPUTS = ["forward", "reverse"]
    os.makedirs("forward", exist_ok=True)
    os.makedirs("reverse", exist_ok=True)
    reads="{dir}/reads.fastq"
else:
    INPUTS = ["forward"]
    os.makedirs("forward", exist_ok=True)
    reads=FASTQ_FILE


rule all:
    input:
        expand("{dir}/output/plots/ScatterReadLengthByAvgEE.png", dir=INPUTS),
        expand("{dir}/output/plots/HeatmapIndexValueByAvgEE.png", dir=INPUTS),
        expand("{dir}/output/plots/HistogramRetainedReadCount.png", dir=INPUTS),
        expand("{dir}/output/SMARTDADA2_InteractiveOutput.html", dir=INPUTS)

if PAIRED:
    rule get_paired_files:
        input:
            fq = FASTQ_FILE
        output:
            forward = INPUTS[0]+"/reads.fastq",
            rev = INPUTS[1]+"/reads.fastq"
        shell:
            "python smartdada2/GetPairedendFiles.py --fq {input.fq} --of_f {output.forward} --of_r {output.rev}"

rule get_subsample:
    input:
        os.path.abspath(reads)
    output:
        "{dir}/subsample.fastq"
    params:
        n = SUBSAMPLE
    shell:
        "python smartdada2/GetSubsamples.py --fq {input} --n {params.n} --of {output}"

rule create_TSVs:
    input:
        "{dir}/subsample.fastq"
    output:
        TrimInfo = "{dir}/output/TrimInfo.tsv",
        SumEEInfo = "{dir}/output/sumEEInfo.tsv"
    params:
        threshold = AVG_Q_SCORE,
        o_mtp = OBV_TRIMMING_MAX,
        a_mtp = MAX_TRIMMING
    shell:
        "python smartdada2/main.py --fq {input} --t_of {output.TrimInfo} --th {params.threshold} --o_mtp {params.o_mtp} --a_mtp {params.a_mtp} --EE_of {output.SumEEInfo}"

rule create_TSVs_checkpoint:
    output:
        touch("{dir}/output/TrimInfo_checkpoint")

rule get_plots:
    input:
        TrimInfo = "{dir}/output/TrimInfo.tsv",
        sumEEInfo = "{dir}/output/sumEEInfo.tsv",
        output_file = "{dir}"
    output:
        "{dir}/output/plots/ScatterReadLengthByAvgEE.png",
        "{dir}/output/plots/HeatmapIndexValueByAvgEE.png",
        "{dir}/output/plots/HistogramRetainedReadCount.png"
    shell:
        "Rscript " + "smartdada2/R_scripts/vizualize.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo} --output_file {input.output_file}"

rule get_interactive_plots:
    input:
        TrimInfo = "{dir}/output/TrimInfo.tsv",
        sumEEInfo = "{dir}/output/sumEEInfo.tsv",
        output_file = "{dir}"
    output:
        "{dir}/output/SMARTDADA2_InteractiveOutput.html"
    shell:
        "Rscript " + "smartdada2/R_scripts/wrapper_RunRender.R --trim_file {input.TrimInfo} --sumEE_file {input.sumEEInfo} --output_file {input.output_file}"

ruleorder: create_TSVs > create_TSVs_checkpoint > get_plots > get_interactive_plots
