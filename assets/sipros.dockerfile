FROM ubuntu:22.04
ENV TZ=Asia/Kolkata \
    DEBIAN_FRONTEND=noninteractive
RUN mkdir /software/

#install python
RUN apt update
RUN apt install python2.7 wget openssl libxml2-dev -y && \
    ln /usr/bin/python2.7 /usr/bin/python

#install python packages
COPY sipros_requirements.txt /software/
RUN cd /software/ && \
    pip install -r sipros_requirements.txt

#install mono
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/ubuntu stable-focal main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt update
RUN apt install mono-devel

#install R
RUN apt install r-base r-base-dev r-cran-tidyverse -y

#intall R packages
COPY install.r /software/
RUN cd /software/ && \
    Rscript install.r &> /software/rsinstall.log

#install sipros
RUN cd /software/ && \
    wget https://github.com/thepanlab/Sipros4/releases/download/4.02/siprosRelease.zip && \
    unzip siprosRelease.zip -d Sipros4 && \
    chmod -R +x Sipros4/


# RUN apt update
# RUN apt install gnupg dirmngr ca-certificates wget aptitude -y
# RUN gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
# RUN echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] http://download.mono-project.com/repo/ubuntu bionic/snapshots/5.18 main" | \
#     tee /etc/apt/sources.list.d/mono-official-stable.list
# RUN apt update
# RUN cd /software/ && \
#     wget https://archive.debian.org/debian/pool/main/libj/libjpeg8/libjpeg8_8d-1+deb7u1_amd64.deb && \
#     apt install ./libjpeg8_8d-1+deb7u1_amd64.deb -y
# RUN aptitude install mono-complete -y

