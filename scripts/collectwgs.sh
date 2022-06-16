#!/bin/bash
ref="/SGRNJ03/randd/user/wuqi/Dna_seq/path/to/project-workdir/resources/genome.fasta"
for bam in `find . -name "*.bam"`
do 
  input="${bam}"
  output="${bam}_collect_wgs_metrics.txt"
  picard CollectWgsMetrics I=$input O=$output R=$ref
done