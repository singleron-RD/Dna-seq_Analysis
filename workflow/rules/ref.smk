"""
rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    script:
        "../../scripts/ref/get_genome.py"
""" 


checkpoint genome_faidx: 
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai"
    log:
        "logs/genome-faidx.log",
    cache: True
    shell:
        "samtools faidx {input} > {output} 2>{log}"


"""
rule genome_dict: 
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "



rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    script:
        "../../scripts/ref/get_known_variation.py"
 


rule remove_iupac_codes: 
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants: 
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    cache: True
    shell:
        "tabix -p vcf {input} 2>{log}"
"""


rule bwa_index: 
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    shell:
        "../../scripts/ref/bwa_index.py"


"""
rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    shell:
        "../../scripts/ref/get_vep_cache.py"



rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    shell:
        "../../scripts/ref/get_vep_plugins.py"
"""