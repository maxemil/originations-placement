From:continuumio/miniconda3:4.5.12
Bootstrap:docker

%labels
    MAINTAINER Max Emil Schön <max-emil.schon@icm.uu.se>
    DESCRIPTION Singularity image containing epa-ng, mafft, raxml-ng and gappa
    VERSION 0.0

%environment
    PATH=/opt/conda/envs/placement/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    apt-get update
    apt-get install -y procps libxtst6
    apt-get clean -y

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a


%test
/opt/conda/envs/placement/bin/mafft --version
/opt/conda/envs/placement/bin/python3 --version
/opt/conda/envs/placement/bin/epa-ng --version
/opt/conda/envs/placement/bin/gappa --help
/opt/conda/envs/placement/bin/trimal --version
