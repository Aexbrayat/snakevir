snakevir

Authors

Antoni Exbrayat (CIRAD) & Etienne Loire (CIRAD) & Serafin Gutierrez (CIRAD)
Purpose:

Metagenomic analysis of viral samples
Steps:

    cleaning
    merging
    filtering
    assembly
    taxonomic annotation

Usage:

    edit config.yaml to precise dataset and dependencies path
    edit snakefile to accomodate read files. Currently set to {sample_name}_1.fastq.gz and {sample_name}_2.fastq.gz
    launch with e.g. :

snakemake -s snakefile --cluster "qsub -q normal.q -V -cwd -pe parallel_smp {threads}" --jobs 100

to execute on a SGE cluster with a maximum of 100 concurrent jobs submitted.
Dependencies
    cutadapt
    bwa
    flash
    megahit
    cap3
    python3
    samtools
    picard tools
    ncbi-blast
    ETE toolkit
