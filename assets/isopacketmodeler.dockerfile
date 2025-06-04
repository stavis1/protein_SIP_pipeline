FROM continuumio/miniconda3:latest

RUN wget https://raw.githubusercontent.com/stavis1/isopacketModeler/refs/heads/main/env/run.yml && conda env create -n isotope_env -f run.yml

