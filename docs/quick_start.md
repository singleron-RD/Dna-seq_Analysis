# Quick start

## Installation
Document can be found [here](installation.md)

## Command-line interface
Assays can be one of:
- ref
- split
- map
- call
- filter
- annotation
- stat
- qc

## Usage Example

1. Make genomeDir
This process can use the ref module to automatically generate the required genomeDir. If necessary, you can also download the required genome file from the bundle at the download address [http://ftp.ensembl.org/pub/](http://ftp.ensembl.org/pub/) ,if you need to download independently, please refer to the example for operation.

### Homo sapiens

```
mkdir homo_sapiens
cd homo_sapiens

curl -L ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | gzip -d > genome.fasta
curl -L ftp://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz | gzip -d > genome.gtf

mkdir -p vep/cache
cd vep
curl -L ftp://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz -o homo_sapiens_vep_108_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_108_GRCh38.tar.gz
```
Note that the above only download genome files, and continue to use the ref module to index genome files.
```
dna-seq ref --config_path ./ --species homo_sapiens --release 108 --build GRCh38 --datatype dna --type all
```


2. Generate scripts

Under your working directory, write a shell script `multi.sh` as

```
conda activate wgs-process
multi_wgs --config_path test_data/\
 --whether_split true\
 --outdir ./\
 --species homo_sapiens\
 --release 108\
 --build GRCh38\
 --thread 8\
 --mod shell\
 --mapfile test_data/file/mapfile
```
`--config_path` Required. The position where the config.yaml is located.Check [multi_wgs.md](./multi_wgs.md) for how to write the config.yaml .

`--whether_split` Required. Whether to split the original sequence according to the index.Please refer to the [index.fasta](../test_data/file/index.fa) file for this index format,all sequences used can be found in [barcode-ADT.fasta](../assets/barcode-ADT.fasta).

`--species` Required. Ensembl species name.

`--release` Required. Ensembl release.

`--build` Required. Genome build.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `run.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/run.sh
```


## Usage

- [multi_wgs.md](./multi_wgs.md)


## Test data
[test_data](../test_data/)