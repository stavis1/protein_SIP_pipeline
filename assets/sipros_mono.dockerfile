FROM ubuntu:18.04
RUN mkdir /software

#install mono 
RUN apt update
RUN apt install ca-certificates gnupg wget unzip -y
RUN gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/ubuntu bionic/snapshots/5.18 main" | tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt update
RUN apt install mono-complete -y

#install sipros
RUN cd /software/ && \
    wget https://github.com/thepanlab/Sipros4/releases/download/4.02/siprosRelease.zip && \
    unzip siprosRelease.zip -d Sipros4 && \
    chmod -R +x Sipros4/
