#!/usr/bin/env bash
set -euxo pipefail
dna-seq_split --config_path ./config;
dna-seq_trimming --config_path ./config;
dna-seq_mapping --config_path ./config;
dna-seq_calling --config_path ./config;
dna-seq_filtering --config_path ./config;
dna-seq_annotation --config_path ./config;
dna-seq_stats --config_path ./config;
dna-seq_qc --config_path ./config