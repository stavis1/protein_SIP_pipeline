FROM python:2.7
RUN mkdir /software

#install sipros
RUN cd /software/ && \
    apt install wget -y && \
    wget https://github.com/thepanlab/Sipros4/releases/download/4.02/siprosRelease.zip && \
    unzip siprosRelease.zip -d Sipros4 && \
    chmod -R +x Sipros4/

#install python packages
COPY sipros_requirements.txt /software/
RUN cd /software/ && \
    pip install -r sipros_requirements.txt