# Start with the barebones image
ARG MINIBASE_IMAGE=gatksv/sv-base-mini:cbb1fc
FROM ${MINIBASE_IMAGE}

ENV HOME=/root


####################################################################################################
# R and all packages
# Steps:
#  1. OS libs
#  2. install base R, see https://cran.r-project.org/bin/linux/ubuntu/README.html#administration-and-maintances-of-r-packages
#  3. install R packages that all/most child images depend on
#  4. remove doc, help, demos and tests (saves ~ 10% space)
#  4. final clean up
COPY install_R_packages.R /opt/
COPY install_bioconductor_packages.R /opt/
ARG R_DEP_TRANSIENT="make gpg-agent gnupg"
ARG R_DEPENDENCIES="software-properties-common build-essential ${R_DEP_TRANSIENT} libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libxml2-dev ca-certificates apt-transport-https openssh-client"
ARG R_RELEASE_VERSION="3.5.1"
ARG SV_BASE_R_PKGS="ggplot2 bbmle emdbook"
ARG SLIM_R_LIB_CMD="find .  -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${R_DEPENDENCIES} && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/" && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 r-base-core="${R_RELEASE_VERSION}"-2bionic \
                 r-base-dev="${R_RELEASE_VERSION}"-2bionic && \
    mkdir -p /tmp/R_pkg_download/ && \
    Rscript -e "update.packages(destdir = \"/tmp/R_pkg_download/\", clean = TRUE, quiet = TRUE)" && \
    cd /opt/ && \
    Rscript --vanilla install_R_packages.R ${SV_BASE_R_PKGS} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy purge \
                 ${R_DEP_TRANSIENT} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

##################################################
# conda and all packages
# Steps:
#  1. OS libs
#  2. install miniconda, enable for all users (instead of current user only), update and cleanup
#  3. install some commonly used conda packages
#  4. final clean up
#  5. special note:  Pysam 0.14.1 - Installing from source prevents the following bug:
#                    python: vcf.c:3482: bcf_update_format: Assertion `!fmt->p_free' failed.
#                    /cromwell_root/script: line 31: 22 Aborted (core dumped) /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py
ARG CONDA_DEP_TRANSIENT="make git wget"
ARG CONDA_DEP="software-properties-common zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev g++ gfortran ${CONDA_DEP_TRANSIENT}"
# versions of bedtools > 2.27.0 seem to have lost the ability to read gzipped files
# pandas 1.0.0 causes problem with bedtools in aggregate.py
ARG PYTHON_PKGS="wheel=0.34.2 bzip2=1.0.8 cython=0.29.14 numpy=1.18.1 pandas=0.25.3 scikit-learn=0.22.1 scipy=1.4.1 intervaltree=3.0.2 matplotlib=3.1.3 natsort=7.0.1 bedtools=2.27.0 pybedtools=0.8.1 pysam=0.16.0.1 bwa=0.7.17 picard=2.23.8 fastq-tools=0.8.3"
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR="/opt/conda"
ARG CONDA_RELEASE="4.6.14"
ARG CONDA_BIN=${CONDA_INSTALL_DIR}/bin
ARG CONDA_CMD=${CONDA_BIN}/conda
ENV PATH=${CONDA_BIN}:$PATH
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${CONDA_DEP} && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR} && \
    ln -s ${CONDA_INSTALL_DIR}/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo "${CONDA_CMD} activate base" >> ${HOME}/.bashrc && \
    /bin/bash -c "source ${HOME}/.bashrc" && \
    conda install -qy \
          -c conda-forge \
          -c bioconda \
          python=3.8 \
          ${PYTHON_PKGS} && \
    git clone https://github.com/pysam-developers/pysam.git /opt/pysam && \
    cd /opt/pysam && \
    git checkout af6ff23322af451b8ab6ea78b53676c45ae0a779 && \
    pip install -e . && \
    rm -rf ${HOME}/.cache/pip && \
    conda clean -tipsy && \
    rm /tmp/miniconda.sh && \
    apt-get -qqy remove --purge ${CONDA_DEP_TRANSIENT} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
####################################################################################################
# Localize github repo containing necessary python and R scripts
# Steps:
#  1. install git
#  2. make directory (/home/git) in which to clone git repo and clone git repo
#  3. copy scripts to /opt/
#  4. remove /home/git
RUN apt -qqy update && \
    apt -qqy install git
RUN mkdir /home/git && \
     cd /home/git && \
     git clone https://github.com/sph17/down-sampling.git
RUN cp /home/git/down-sampling/scripts/* /opt/
RUN rm -rf /home/git/

