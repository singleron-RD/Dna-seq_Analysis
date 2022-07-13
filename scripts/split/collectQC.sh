#!/bin/bash
ref="../../resources/genome.fasta"
refflat="../../resources/genome.refFlat"
for bam in `find . -name "*.bam"`
do 
  input="${bam}"
  output_picard_wgs="${bam}_collect_wgs_metrics.txt"
  output_samtools="${bam}_mapping.txt"
  picard CollectWgsMetrics I=$input O=$output_picard R=$ref;samtools stats $input > $output_samtools;\
  python ../../scripts/split/collect_feature.py --bam $input --outdir ./ --refflat $refflat --sample $input
done