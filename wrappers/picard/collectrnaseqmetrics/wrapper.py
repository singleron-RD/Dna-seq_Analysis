__author__ = "wuqi"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

region_log = snakemake.output

shell(
    "picard CollectRnaSeqMetrics "  # Tool and its subcommand
    "-Xmx20G"
    "-XX:ParallelGCThreads=8"
    "I={snakemake.input.bam} "  # Input bam(s)
    "O={region_log} "  # Output bam
    "REF_FLAT {snakemake.input.refflat}"
    "STRAND=NONE"
    "VALIDATION_STRINGENCY=SILENT"
    "{log}"  # Logging
)