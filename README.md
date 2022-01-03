snakevir

Authors

Antoni Exbrayat (CIRAD) & Etienne Loire (CIRAD) & Serafin Gutierrez (CIRAD)

Purpose:
Metagenomic analysis of viral shotgun NGS samples. This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/). 

## Dependencies  
  - bioawk
  - biopython
  - blast
  - bwa
  - cap3
  - csvkit
  - cutadapt
  - diamond
  - entrez-direct
  - ete3
  - flash
  - megahit
  - pandas
  - picard
  - python
  - r-base
  - samtools
  - seqtk
  - snakemake

The conda environment manager can be used to install python , snakemake and all the required tools and dependencies into a single environment in a way such that reproducibility is ensured. 

Note: Conda must be installed on the system. For help with setting up conda, please see [miniconda](https://docs.conda.io/en/latest/miniconda.html).

To create and activate the conda environment with the environment.yml provided , use :
```
conda env create -f environment.yml
conda activate snakevir
```

## Steps:

    cleaning
    merging
    filtering
    assembly
    taxonomic annotation

## Usage:
    edit config.yaml to precise dataset and dependencies path
    edit snakefile to accomodate read files. Currently set to {sample_name}_1.fastq.gz and {sample_name}_2.fastq.gz
    launch with e.g. :

snakemake -s snakefile --cluster "qsub -q normal.q -V -cwd -pe parallel_smp {threads}" --jobs 100

to execute on a SGE cluster with a maximum of 100 concurrent jobs submitted.
