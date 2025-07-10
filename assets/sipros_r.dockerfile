FROM bioconductor/bioconductor_docker:3.21-R-4.5.1
RUN mkdir /software

#install sipros
RUN cd /software/ && \
    apt install wget -y && \
    wget https://github.com/thepanlab/Sipros4/releases/download/4.02/siprosRelease.zip && \
    unzip siprosRelease.zip -d Sipros4 && \
    chmod -R +x Sipros4/
