# 使用手册

## 需要拆分的样本

如果需要根据index去拆分原始序列，首先需要写入config/mafile文件；
随后将config/config.yaml中split_fq设置为split_fq: true；
最后直接使用 snakemake 流程进行常规分析。

所需mapfile文件：

```shell
$cat config/mapfile
sample_name	fq_outdir	fa	fq1	fq2	mod
test1	./test_data/fastq/	./test_data/file/index1.fa	./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq	fq2 
test2	./test_data/fastq/	./test_data/file/index2.fa	./test_data/fastq/test_R1.fastq	./test_data/fastq/test_R2.fastq	fq2 
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

mod为拆分使用的fq序列。



## 不需要根据index拆分的样本
首先编辑好 config/samples.tsv 及 config/units.tsv 文件。
随后将config/config.yaml中split_fq设置为split_fq: false；
最后直接使用 snakemake 流程进行常规分析。


## snakemake 使用说明
修改config/config.yaml文件中的前两行，保证所用的 samples.tsv 和 units.tsv 与将要进行分析的样本一致。
所使用的参考基因组保存在 resources 目录下，可提前下载好，如未提前下载，该流程会自动下载部署，所需时间较长，并且由于网络原因可能存在下载不完全的问题，此点需要注意。
此流程适用于 release-98_homo_sapiens_GRCh38 参考基因组，如需更换参考数据库及物种，请修改 config/config.yaml 中的ref信息。

运行前先检查流程是否能正常运行：

```shell
 snakemake -s workflow/Snakefile --cores 10 -np 
```

本地直接运行： 

```shell
snakemake -s workflow/Snakefile --cores 10 
#--cores 指定所用的cpu数
```


提交至集群运行： 

```shell
snakemake -s workflow/Snakefile  --cluster "qsub -cwd -V -l vf=10g,p=3 -q randd.q -o o.logs -e e.logs" --jobs 1 
#jobs 指定能同时提交的任务数
```
### snakemake 结果说明

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



## 注意事项

### 创建conda环境

流程运行前，请确保环境中已安装相关分析软件

```shell
conda env create -f environment.yaml;pip install requirements.txt
```

### 目录一致性

此流程的目录结构请保证如下所示：
.
├── config
├── js
├── resources
├── scripts
│   ├── calling
│   ├── filtering
│   ├── mapping
│   ├── qc
│   ├── ref
│   ├── split
│   └── stats
└── workflow

在当前目录下运行 `snakemake -s workflow/Snakefile --cores 10`后，会生成 ./results 目录以及 ./.snakemake 目录。
此后需要运行其他样本请修改 ./results 目录名称，防止新生成的结果覆盖原结果文件。

