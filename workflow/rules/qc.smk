rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html=temp("results/qc/fastqc/{sample}-{unit}.html"),
        zip=temp("results/qc/fastqc/{sample}-{unit}.zip"),
    log:
        "logs/fastqc/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.fastqc.benchmark.txt",
    script:
        "../../scripts/qc/fastqc.py"


rule samtools_stats:
    input:
        "results/recal/{sample}-{unit}.bam",
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt",
    log:
        "logs/samtools-stats/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.samtools_stats.benchmark.txt",
    shell:
        "samtools stats {input} > {output} 2>{log}"


rule multiqc:
    input:
        expand(
            [
                "results/qc/samtools-stats/{u.sample}-{u.unit}.txt",
                "results/qc/fastqc/{u.sample}-{u.unit}.zip",
                "results/qc/dedup/{u.sample}-{u.unit}.metrics.txt",
            ],
            u=units.itertuples(),
        ),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    script:
        "../../scripts/qc/multiqc.py"