FROM continuumio/miniconda3:latest

RUN conda create -n sipros_env -c conda-forge -c bioconda sipros


