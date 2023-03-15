## Usage
- The original sequence needs to be split according to the index:
```
    multi_wgs --config_path test_data/\
              --whether_split true\
              --mapfile test_data/file/mapfile\
              --outdir ./\
              --species homo_sapiens\
              --release 108\
              --build GRCh38\
              --thread 8
              
```
Note: Please set the `whether_split` parameter to `true` and write the `mapfile` parameter.

- No need to split the original sequence according to the index:
```
    multi_wgs --config_path test_data/\
              --whether_split false\
              --outdir ./\
              --species homo_sapiens\
              --release 108\
              --build GRCh38\
              --thread 8
              
```
Note: Please set the `whether_split` parameter to `false`.

Note: If use the same reference genome, only need to run the ref module once, and then can adjust the `steps_not_run` parameter to `ref`.
```
    multi_wgs --config_path test_data/\
              --whether_split true\
              --mapfile test_data/file/mapfile\
              --outdir ./\
              --thread 8\
              --steps_not_run ref

```

## Features
### ref
- Create a genome reference directory.

### split

- Split the original sequence data into sequence data of single cell according to index.

### trim
- Remove adapter sequence and perform quality control in R1 and R2 reads with trimmomatic. Default parameters includes:
    - LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### map
- Align paired-end reads to the reference genome with bwa mem.

- Recalibrate the base mass with gatk.

### call
- SNP Calling with gatk.

### filter
- Filter mutation files.

### annotation
- Variation note with VEP.

### stat
- Statistics of some mutation information.

### qc
- Make statistics on the quality of the original sequence, the comparison repetition rate and other information with multiqc.


## Output files
### ref
- genome.fasta

### split
- The split sequence according to the index.
    - `{sample_1}_{index_1}_R1.fastq`,`{sample_1}_{index_1}_R2.fastq` Paired-end read sequences of sample1 split according to index1.

    - `{sample_1}_not_R1.fastq`,`{sample_1}_not_R2.fastq` Invalid sequence of sample1 split according to index.

### trim
- `01.trimmed/{sample-1}.1.fastq.gz`,`01.trimmed/{sample-1}.2.fastq.gz` paired fastq of sample1.

- `01.trimmed/{sample-1}.1.unpaired.fastq.gz`,`01.trimmed/{sample-1}.2.unpaired.fastq.gz` unpaired fastq of sample1.

- `01.trimmed/{sample-1}.trimlog.txt` log of trimmomatic.

### map
- `02.mapped/{sample}-{1}.sorted.bam` Sorted BAM file contains Uniquely Mapped Reads.

- `02.mapped/{sample}-{1}_chrom.coverage.txt` Coverage and depth information from sequence alignment of sample to genome.

- `02.mapped/{sample}-{1}_mapping.txt` Bam information file of samtools stats statistics.

- `02.mapped/merge.tsv` The information file of all samples is consolidated.

- `03.dedup/{sample}-{1}.bam` BAM file after remove duplicates reads. 

- `04.recal/{sample}-{1}.bam` Bam file after base mass correction.

### call
- `06.genotyped/all.vcf.gz` VCF file after mutation detection.

### filter
- `07.filted/all.indels.vcf.gz` VCF file after extract indel.

- `07.filted/all.snvs.vcf.gz` VCF file after extract snp.

- `07.filted/all.vcf.gz` VCF file after merging snp vcf and indel vcf.

### annotation
- `08.annotated/all.vcf.gz` VCF file after vep annotation.

- `09.stats/all_sgr.html` Statistics of all samples.Include General statistics,Variant classes,Consequences (most severe),Variants by chromosome and other information.

- `10.tables/calls.tsv.gz` Extract DP, AD and other information of all samples.

- `11.split/{sample}/{sample}*` Information of individual samples.

### stat
- `12.plots/allele-freqs.svg` Scatter chart of mutation frequency of all samples.

- `12.plots/depths.svg` Scatter chart of reads depth of all samples.

### qc
- `qc/multiqc.html` Aggregate results from bioinformatics analyses across many samples into a single report.

### logs
Log file directory of process output.


## Arguments
`--config_path` The position where the config.yaml is located. 
The format of the config.yaml file is as follows.
```
$cat config.yaml
outdir: ./results 
genomedir: ./homo_sapiens
```
Note that if `whether_split` is true, the required units.tsv file will be automatically generated,if `whether_split` is false, you need to write the units.tsv file into the config_path directory.The format of the units.tsv file is as follows.
```
$cat units.tsv
sample	unit	platform	fq1	fq2
test1_index_1	1	ILLUMINA	test_data/split_fq/test1_index_1_R1.fastq	test_data/split_fq/test1_index_1_R2.fastq
test1_index_2	1	ILLUMINA	test_data/split_fq/test1_index_2_R1.fastq	test_data/split_fq/test1_index_2_R2.fastq
test2_index_3	1	ILLUMINA	test_data/split_fq/test2_index_3_R1.fastq	test_data/split_fq/test2_index_3_R2.fastq
test2_index_4	1	ILLUMINA	test_data/split_fq/test2_index_4_R1.fastq	test_data/split_fq/test2_index_4_R2.fastq
```

`--mod` Which type of script to generate, `sjm` or `shell`.

`--queue` Only works if the `--mod` selects `sjm`.

`--rm_files` Remove redundant bam files after running.

`--steps_run` Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `split` and `trim`, 
use `--steps_run split,trim`.

`--steps_not_run` Steps not to run. Multiple Steps are separated by comma. For example, if you want  not to run `ref` and `split`, 
use `--steps_not_run ref,split`.

`--whether_split` Whether to split the original sequence according to the index.

`--outdir` The directory of script output.

`--thread` Thread to use.

`memory` The memory size used, the default is GB.

`--species` Ensembl species name.

`--release` Ensembl release.

`--build` Genome build.
 
`--datatype` Sequence types.

`--chromosome` Select a specific chromosome for analysis.

`--type` Ensembl VCF (Variant Call Format) files types.

`--index_prefix` The prefix of the generated index files.(same as fasta name)

`--algorithm` BWT construction algorithm: bwtsw, is or rb2 [auto]

`--mapfile` The necessary mapfile file location for splitting the original sequence.
The format of the mapfile file is as follows.
```
$cat mapfile
sample_name	fq_outdir	fa	fq1	fq2
test1	test_data/split_fq	test_data/file/index1.fa	test_data/fastq/test_R1.fastq	test_data/fastq/test_R2.fastq
test2	test_data/split_fq	test_data/file/index2.fa	test_data/fastq/test_R1.fastq	test_data/fastq/test_R2.fastq
```
Note that this parameter only takes effect when the `whether_split` parameter is true, which is required to run the split step.
`--trim_param` Additional parameters for the called software. Need to be enclosed in quotation marks.For example, `--{software}_param 'param1:value1 param2:value2'`,`--trim_param 'LEADING:3 TRAILING:3'`.

`--remove_duplicates` Delete the duplicate sequence after comparison.

`--intervals` One or more genomic intervals over which to operate.This argument may be specified 0 or more times.

`--interval_padding` Amount of padding (in bp) to add to each interval you are including.

`--rm_files` Remove redundant bam files after running.

`--filtertype` Types of filter.

`--vep_param` Additional parameters for the called software. Need to be enclosed in quotation marks.For example, `--{software}_param '--param1 value1 --param2 value2'`,`--vep_param '--sift b'`.

`--vep_plugins_param` Add any plugin from [https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html).Use named plugin.Multiple plugins can be used by supplying the --vep_plugins_param flag multiple times. For example, `--vep_plugins_param 'LoFtool'`.
