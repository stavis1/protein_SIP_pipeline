FROM python:2.7

RUN mkdir /software/

#install mono
RUN apt update
RUN apt install gnupg dirmngr ca-certificates wget -y
RUN gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/debian stable-buster main" | \
    tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt update
RUN apt install mono-complete -y

#install python packages
COPY sipros_requirements.txt /software/
RUN cd /software/ && \
    pip install -r sipros_requirements.txt

#install R
RUN gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null
RUN apt update
RUN apt install r-base r-base-dev -y

#intall R packages
COPY install.r /software/
RUN cd /software/ && \
    Rscript install.r

#install sipros
RUN cd /software/ && \
    git clone https://github.com/thepanlab/Sipros4 && \
    chmod -R +x Sipros4/


