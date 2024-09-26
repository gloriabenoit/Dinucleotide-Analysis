# Dinucleotide frequencies analysis in reference genomes
May 2023 (L3 BI)

## Introduction

This project aims to compute dinucleotide frequencies in several reference genomes and compare them with theoretical values.
This is inspired by a blog post written by Guillaume Devailly titled [Fréquences des dinucléotides dans le génome d'organismes modèles](https://bioinfo-fr.net/frequences-des-dinucleotides-dans-le-genome-dorganismes-modeles).

## Setup

To install the algorithm and its dependencies, you need to perform the following steps:

### Clone the repository

```bash
git clone https://github.com/gloriabenoit/Dinucleotide-Analysis.git

cd Dinucleotide-Analysis
```

### Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Create a Conda environment

```bash
conda env create -f environment.yml
```

### Activate the Conda environment

```bash
conda activate dinucl-freq
```

## Usage (command line interface)

```bash
python src/GloriaBENOIT_DinucleotideAnalysis.py fasta_file_1 ... fasta_file_n
```

## Data used

Since reference genomes are quite large, I have not included them. However you can find the download links in `doc/all_genomes.txt`.
