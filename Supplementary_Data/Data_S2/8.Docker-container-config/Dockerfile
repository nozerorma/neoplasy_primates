# Dockerfile for building immutable environment
# used for Primate Neoplasy Masters Project

FROM rocker/r-ver:4.3.1

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

LABEL image.author.name "Miguel Ramon"
LABEL image.author.email "miguel.ramon@upf.edu"

ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2023.09.0+463
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default

# Remove ulimit set by Rstudio dist
#RUN ulimit -n 1048576 # to build in podman remove this

# Download some basic debian tools + procps (Nextflow dependency).
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -qq -y --no-install-recommends \
        software-properties-common \
        dirmngr \
        ed \
        gpg-agent \
		less \
		locales \
		vim \
		ca-certificates \
        libglpk-dev \
        libxt6 \
        python3 \
        pip \
        curl \
        ssh \
        tree
        
# Install RStudio stuff
RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh
RUN /rocker_scripts/install_tidyverse.sh

# Set WD and create build directory structure.
RUN mkdir -p /home/rstudio/NEOPLASY_PRIMATES
WORKDIR /home/rstudio/NEOPLASY_PRIMATES

# Install some stuff for clustering
RUN pip3 install \
        treecluster \
        treeswift \
        niemads

EXPOSE 8787

CMD ["/init"]