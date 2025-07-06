FROM python:3.11.9

RUN mkdir /software/
COPY isopacketmodeler_requirements.txt /software/
RUN cd /software/ && \
    pip install -r isopacketmodeler_requirements.txt && \
    git clone https://github.com/stavis1/isopacketModeler && \
    cd isopacketModeler && \
    pip install ./

RUN mkdir /scripts/
COPY sipros2IPM.py /scripts/
COPY split_peptides.py /scripts/
COPY merge_peptides.py /scripts/
COPY formula_parser.py /scripts/