rule get_genome:
    output:
        "resources/homo_sapiens/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "file:wrappers/reference/ensembl-sequence"



checkpoint genome_faidx: 
    input:
        "resources/homo_sapiens/genome.fasta",
    output:
        "resources/homo_sapiens/genome.fasta.fai"
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "file:wrappers/samtools/faidx"


rule genome_dict: 
    input:
        "resources/homo_sapiens/genome.fasta",
    output:
        "resources/homo_sapiens/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/homo_sapiens/genome.fasta.fai",
    output:
        vcf="resources/homo_sapiens/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "file:wrappers/reference/ensembl-variation"


rule remove_iupac_codes: 
    input:
        "resources/homo_sapiens/variation.vcf.gz",
    output:
        "resources/homo_sapiens/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants: 
    input:
        "resources/homo_sapiens/variation.noiupac.vcf.gz",
    output:
        "resources/homo_sapiens/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "file:wrappers/tabix"


rule bwa_index: 
    input:
        "resources/homo_sapiens/genome.fasta",
    output:
        multiext("resources/homo_sapiens/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "file:wrappers/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/homo_sapiens/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "file:wrappers/vep/cache"



rule get_vep_plugins:
    output:
        directory("resources/homo_sapiens/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "file:wrappers/vep/plugins"
