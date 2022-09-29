# 使用手册

## 软件安装

```shell
conda env create -f environment.yaml;pip install requirements.txt;python setup.py
```

## 需要拆分的样本

如果需要根据index去拆分原始序列，首先需要写入config/mafile文件;
随后运行 `dna-seq_split --config_path ./config`;
再运行 `run_dna.sh` 脚本进行常规分析。

```shell
$ cat run_dna.sh
#!/usr/bin/env bash
set -euxo pipefail

dna-seq_split --config_path ./config;   # need split
dna-seq_trimming --config_path ./config;
dna-seq_mapping --config_path ./config;
dna-seq_calling --config_path ./config;
dna-seq_filtering --config_path ./config;
dna-seq_annotation --config_path ./config;
dna-seq_stats --config_path ./config;
dna-seq_qc --config_path ./config
```


所需mapfile文件：

```shell
$cat config/mapfile
sample_name	fq_outdir	fa	fq1	fq2
test1	./test_data/fastq/	./test_data/file/index1.fa	./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq
test2	./test_data/fastq/	./test_data/file/index2.fa	./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq
test3	./test_data/fastq/	./test_data/file/index.fa	./test_data/fastq/test_R1.fastq,./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq,./test_data/fastq/test_R2.fastq
```

其中fa为使用的 index 文件

```shell
$cat ./test_data/file/index1.fa
 >index1
 AACCGCGGT
 >index2
 GGTTATAAT
 >index3
 CCAAGTCCT
 >index4
 TTGGACTTT
```


## 不需要根据index拆分的样本
首先编辑好 config/units.tsv 文件;
再直接运行 `run_dna.sh` 脚本进行常规分析。


## config.yaml 文件说明
必须写入项:units, outdir, resource。拆分写入项:mapfile。
输出路径：outdir ;
使用参考基因组位置:resource。
所使用的参考基因组保存在 resources 目录下，使用 dna-seq_ref --config_path ./config 获取所需参考基因组信息。
此流程适用于 release-98_homo_sapiens_GRCh38 参考基因组，如需更换参考数据库及物种，请修改 config/config.yaml 中的ref信息。


###  结果说明

部分结果说明：
`trimmed` 将Illumina平台FASTQ序列中的Adapter去除，并根据基础质量值对FASTQ序列文件进行裁剪后的序列文件。
`recal` 删除重复序列并进行碱基质量重校正后的Bam文件。
`genotyped` 经过GATK分析后的文件。
`filtered` 进行过滤后的vcf文件。
`annotation` 使用Vep注释后的vcf文件。
`stats` 指标汇总和统计显示。
`qc` 一些质控文件。
`tables` 从注释后的VCF中提取特定指标后的vcf文件。
`split` 将个样本拆分出来后的文件夹。


## 书写规范

config/config.yaml,mapfile,units.tsv书写时保证全部字符为英文，分隔符使用规范，以防识别不出报错。路径尽量设置为绝对路径。