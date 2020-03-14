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
        python3-dev \
        wget \
        git \
        curl \
        make \
        unzip \
        openjdk-11-jdk \
        samtools \
        genometools \
        python3-numpy \
        python3-matplotlib \
        python3-pysam \
        python3-htseq \
        python3-pip \
    && apt-get clean

#install hisat2
RUN curl ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -o hisat2-2.1.0.zip \
    && unzip hisat2-2.1.0.zip -d /usr/local/bin/ \
    && rm hisat2-2.1.0.zip

#add hisat2 dir to path
ENV PATH $PATH:/usr/local/bin/hisat2-2.1.0/

#install trimmomatic
RUN curl http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -o trimmomatic-0.39.zip \
    && unzip trimmomatic-0.39.zip -d /usr/local/bin/ \
    && rm trimmomatic-0.39.zip

#add trimmomatic dir to path
ENV PATH $PATH:/usr/local/bin/Trimmomatic-0.39/

#install gffread
RUN git clone https://github.com/gpertea/gclib /usr/local/bin/gclib \
    && git clone https://github.com/gpertea/gffread /usr/local/bin/gffread \
    && cd /usr/local/bin/gffread/ \
    && make release

#add trimmomatic dir to path
ENV PATH $PATH:/usr/local/bin/gffread/

#change working directory
WORKDIR /opt/

#chmod of /opt/
RUN chmod -R 777 /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]