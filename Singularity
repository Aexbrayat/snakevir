Bootstrap: docker
From: continuumio/miniconda3

# sudo singularity build snakemake Singularity

%files
    snakemake_tutorial.scif
    Snakefile
    config.yaml
    complet_taxo.sh

%environment
    PATH=/opt/conda/bin:$PATH
    export PATH

%post
    apt-get update && apt-get -y install build-essential valgrind time python-numpy python-qt4 python-lxml python-six python-dev
    /opt/conda/bin/conda config --add channels defaults
    /opt/conda/bin/conda config --add channels conda-forge
    /opt/conda/bin/conda config --add channels bioconda

    # Install scif and scif-apps
    /opt/conda/bin/pip install scif
    /opt/conda/bin/scif install /snakemake_tutorial.scif

    # Install snakemake
    /opt/conda/bin/pip install snakemake==4.4.0
    /opt/conda/bin/pip install docutils==0.14
    /opt/conda/bin/pip install biopython
    /opt/conda/bin/pip install pandas

%runscript
    PATH=/opt/conda/bin:$PATH
    export PATH
    exec scif "$@"
