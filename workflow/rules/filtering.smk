rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/all.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.vcf.gz"),
    params:
        extra=get_vartype_arg,
    log:
        "logs/gatk/selectvariants/{vartype}.log",
    benchmark:
        "results/benchmarks/{vartype}.select_calls.benchmark.txt",
    script:
        "../../scripts/filtering/select_calls.py"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.hardfiltered.vcf.gz"),
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    benchmark:
        "results/benchmarks/{vartype}.hard_filter_calls.benchmark.txt",
    script:
        "../../scripts/filtering/hard_filter_calls.py"


rule recalibrate_calls:
    input:
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf="results/filtered/all.{vartype}.recalibrated.vcf.gz",
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"],
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log",
    benchmark:
        "results/benchmarks/{vartype}.recalibrate_calls.benchmark.txt",
    script:
        "../../scripts/filtering/recalibrate_calls.py"


rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    benchmark:
        "results/benchmarks/merge-filtered.merge_calls.benchmark.txt",
    script:
        "../../scripts/calling/merge_variants.py"
