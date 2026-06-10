FROM python:3.11.15-trixie

RUN wget https://github.com/percolator/percolator/releases/download/rel-3-09/percolator-v3-09-linux-amd64.deb && \
	apt-get update -y && \
	apt-get upgrade -y && \
	apt-get install -y libboost-filesystem1.83.0 && \
	apt-get install -y ./percolator-v3-09-linux-amd64.deb

RUN mkdir /software/
COPY percolator_requirements.txt /software/
RUN pip install -r /software/percolator_requirements.txt
COPY percolator.py /software/

