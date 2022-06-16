rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.trim_reads_se.benchmark.txt",
    wrapper:
        "file:wrappers/trimmomatic/se"


rule trim_reads_pe: 
    input:
        unpack(get_fastq),
    output:
        r1=("results/trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="results/trimmed/{sample}-{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.trim_reads_pe.benchmark.txt",
    wrapper:
        "file:wrappers/trimmomatic/pe"


rule map_reads: 
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.map_reads.benchmark.txt",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "file:wrappers/bwa/mem"


rule mark_duplicates: 
    input:
        "results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=temp("results/dedup/{sample}-{unit}.bam"),
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.mark_duplicates.benchmark.txt",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "file:wrappers/picard/markduplicates"



rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.fasta",
        dict="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.dict",
        known="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/variation.noiupac.vcf.gz",
        known_idx="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.recalibrate_base_qualities.benchmark.txt",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    resources:
        mem_mb=4096, #old:mem_mb=1024,
    wrapper:
        "file:wrappers/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.fasta",
        dict="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam=protected("results/recal/{sample}-{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    benchmark:
        "results/benchmarks/{sample}-{unit}.apply_base_quality_recalibration.benchmark.txt",
    params:
        extra=get_regions_param(),
    resources:
        mem_mb=4096, #old:mem_mb=1024,
    wrapper:
        "file:wrappers/gatk/applybqsr"

rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    benchmark:
        "results/benchmarks/{prefix}.samtools_index.benchmark.txt",
    wrapper:
        "file:wrappers/samtools/index"