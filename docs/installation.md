# Software Requirements
conda
git
# Installation
- Clone repo
```
git clone https://github.com/singleron-RD/Dna-seq_Analysis.git
```
- Create conda environment and install conda packages. It is recommended to use mamba (which is a faster replacement for Conda)
```
cd Dna-seq_Analysis
conda install mamba
mamba create -n wgs-process -y --file conda_pkgs.txt
```
- Install wgs-process
Make sure you have activated the wgs-process conda environment before running python setup.py.
```
conda activate wgs-process
pip install -r requirements.txt
python setup.py develop
```