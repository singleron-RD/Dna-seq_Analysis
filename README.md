# Usage
`Full operation` snakemake -s workflow/Snakefile --use-conda --cores 10 
`Pre-operation inspection` snakemake -s workflow/Snakefile --use-conda --cores 10 -np
`The cluster running` snakemake -s workflow/Snakefile --cluster "qsub -cwd -V -l vf=30g,p=6 -q randd.q -o o.logs -e e.logs" --jobs 6


# Notice
* Your environment must already have Snakemake installed.
* This process is currently used to compare sequencing data from the human reference genome.Reference genome link at `http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/` .


# Description of the snakefile file
`workdir` The working directory for this run. This is where the result files are generated.You need to modify it after running, otherwise the generated results directory will be overwritten.Nothing else needs to be changed.
`configfile and include` If you want to copy it to your own work environment, just change the absolute path to this configuration.


# results
`trimmed` The Adapter in the FASTQ sequence of Illumina platform was removed, and the FASTQ sequence file was trimmed according to the base quality value.
`recal` Bam files with duplicates sequences removed and base quality recalibrated.
`genotyped` VCF file after GATK calling SNP.
`filtered` Filtered VCF files.
`annotation` Vep annotated file.
`stats` Indicators summary and statistics display.
`qc` Quality control documents.
`tables` Extract files for specific metrics from annotated VCF.


## Split fastq
# Usage
python split_fastq.py --mapfile test_data/file/mapfile;
snakemake -s workflow/Snakefile --use-conda --cores 10

## Split vcf
# Usage
python split_vcf.py --qc no


#### 使用手册
### 分步运行
## 需要拆分的样本
如果需要根据index去拆分原始序列，那么需要使用 split_scripts/split_fastq.py 将原始序列拆分为含有各 index 的序列。
# split_fastq.py 脚本使用说明
运行脚本： python split_fastq.py --mapfile test_data/file/mapfile
所需mapfile文件：
`$cat mapfile`
sample_name	fq_outdir	config_outdir	fa	fq1	fq2	mod
test	./test_data/fastq/	./config/	./test_data/file/index.fa	./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq	fq2 
其中fa为使用的 index 文件
`$cat ./test_data/file/index.fa`
 ·>index1
 ·AACCGCGGT
 ·>index2
 ·GGTTATAAT
 ·>index3
 ·CCAAGTCCT
 ·>index4
 ·TTGGACTTT
mod为拆分使用的fq序列。

随后使用 snakemake 流程进行常规分析。
# snakemake 使用说明
修改config/config.yaml文件中的前两行，保证所用的 samples.tsv 和 units.tsv 与将要进行分析的样本一致。

运行前先检查流程是否能正常运行： snakemake -s workflow/Snakefile --use-conda --cores 10 -np 。
本地直接运行： snakemake -s workflow/Snakefile --use-conda --cores 10 (--cores 指定所用的cpu数) 。
提交至集群运行： snakemake -s workflow/Snakefile  --use-conda --cluster "qsub -cwd -V -l vf=30g,p=15 -q randd.q -o o.logs -e e.logs" --jobs 15 (jobs 指定能同时提交的任务数)。


## 拆分报告
此步需在运行完 snakemake 常规分析后运行。
如果是多样本进行的分析，报告将会整合成一份，需要各自的报告需要另外运行 split_scripts/split_vcf.py 。此脚本可以将多样本整合的报告以及vcf文件拆分开。
运行此脚本可以选择是否输出质控文件，质控文件在results/recal目录下。
# split_vcf.py 脚本使用说明
python split_vcf.py --qc no


### 全部写入sh运行
`$cat run.sh`
#!bin/bash
python split_fastq.py --mapfile test_data/file/mapfile;
snakemake -s workflow/Snakefile --use-conda --cores 10;
python split_vcf.py --qc no