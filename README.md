snakevir

Authors

Antoni Exbrayat (CIRAD) & Etienne Loire (CIRAD) & Serafin Gutierrez (CIRAD)

Purpose:
Metagenomic analysis of viral shotgun NGS samples. This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/). 
## Step:
  - Cleaning
  - Merging
  - Filtering
  - De novo sequence assembly
  - Mapping
  - Homology search protein databases
  - Homology search nucleotide databases
  - Taxonomic annotation
  - Taxonomy refining
  - Viral hosts search    
  
## Dependencies  
```
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
``` 
The conda environment manager can be used to install python , snakemake and all the required tools and dependencies into a single environment in a way such that reproducibility is ensured. 

Note: Conda must be installed on the system. For help with setting up conda, please see [miniconda](https://docs.conda.io/en/latest/miniconda.html).

To create and activate the conda environment with the environment.yml provided , use :
```
conda env create -f environment.yml
conda activate snakevir
```


## Usage:
  Snakemake supports a separate configuration file for execution on a cluster. A cluster config file  cluster.json is provided , it allows you to specify cluster   submission parameters outside the Snakefile. The cluster config is contains all parameters with match names of rules in the Snakefile.
  
  edit config.yaml to precise dataset and dependencies path, accomodate read files names , threads allocated to the rules (according to cluster.json).
  
  launch with e.g. :

    snakemake -s snakefile  -j 100  --cluster-config cluster.json --cluster "sbatch -p {cluster.queue} -N {cluster.queue} -c {cluster.cpu_task} --mem {cluster.mem} -e {cluster.error} -o {cluster.log} "  --printshellcmd --rerun-incomplete  --reason --dryrun

  to execute on a SLURM cluster with a maximum of 100 concurrent jobs submitted, eventually modify the command accordingly with your job scheduler.
