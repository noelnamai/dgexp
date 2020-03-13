#specifying the base image
FROM ubuntu:19.04

#add metadata to image e.g. name and email of the person who maintains the file
LABEL Maintainer="Noel Namai"
LABEL Email="noelnamai@yahoo.com"

#set system arguments and enviromental variables
ARG DEBIAN_FRONTEND=noninteractive

#http://bugs.python.org/issue19846: setting "LANG=C" on a Linux system *fundamentally breaks Python 3*
ENV LANG C.UTF-8

#update ubuntu and install other necessary software
RUN apt-get -y update --fix-missing \
    && apt-get install -y --allow-unauthenticated \
        software-properties-common \
        build-essential \
        wget \
        curl \
        unzip \
        samtools \
        genometools \
        python-numpy \
        python-matplotlib \
        python-pysam \
        python-htseq \
        python-pip \
    && apt-get clean

#install hisat2
RUN curl ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -o hisat2-2.1.0.zip \
    && unzip hisat2-2.1.0.zip -d /usr/local/bin/ \
    && rm hisat2-2.1.0.zip

#add hisat2 dir to path
ENV PATH $PATH:/usr/local/bin/hisat2-2.1.0/

#change working directory
WORKDIR /opt/

#chmod of /opt/
RUN chmod -R 777 /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]