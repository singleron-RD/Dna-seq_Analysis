if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.fasta",
        idx="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.dict",
        known="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/variation.noiupac.vcf.gz",
        tbi="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/variation.noiupac.vcf.gz.tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("results/called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    benchmark:
        "results/benchmarks/{sample}.{contig}.call_variants.benchmark.txt",
    params:
        extra=get_call_variants_params,
    wrapper:
        "file:wrappers/gatk/haplotypecaller_v59"


rule combine_calls:
    input:
        ref="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.fasta",
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf="results/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    benchmark:
        "results/benchmarks/{contig}.combine_calls.benchmark.txt",
    wrapper:
        "file:wrappers/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="/SGRNJ06/randd/public/wgs_ref/homo_sapiens/genome.fasta",
        gvcf="results/called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    benchmark:
        "results/benchmarks/{contig}.genotype_variants.benchmark.txt",
    wrapper:
        "file:wrappers/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    benchmark:
        "results/benchmarks/merge-genotyped.merge_variants.benchmark.txt",
    wrapper:
        "file:wrappers/picard/mergevcfs"
