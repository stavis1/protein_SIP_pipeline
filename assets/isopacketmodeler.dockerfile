FROM python:3.11.9

RUN mkdir /software/
COPY isopacketmodeler_requirements.txt /software/
RUN pip install -r /software/isopacketmodeler_requirements.txt

RUN mkdir /scripts/
COPY sipros2IPM.py /scripts/
COPY split_peptides.py /scripts/
COPY merge_peptides.py /scripts/
COPY formula_parser.py /scripts/