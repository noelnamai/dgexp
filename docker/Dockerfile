#specifying the base image
FROM ubuntu:19.04

#add metadata to image e.g. name and email of the person who maintains the file
LABEL Maintainer="Noel Namai"
LABEL Email="noelnamai@yahoo.com"

#set system arguments and enviromental variables
ARG DEBIAN_FRONTEND=noninteractive

#setting "LANG=C" on linux fundamentally breaks python 3
#http://bugs.python.org/issue19846
ENV LANG C.UTF-8

#update ubuntu and install other necessary software
RUN apt-get -y update \
    && apt-get install -y software-properties-common \
    && apt-get -y update \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/" \
    && add-apt-repository ppa:openjdk-r/ppa \
    && apt-get -y update --fix-missing \
    && apt-get install -y --allow-unauthenticated \
        build-essential \
        libssh2-1-dev \
		libxml2-dev \
		libcairo2-dev \
		libxt-dev \
		libssl-dev \
		libv8-dev \
        libmagick++-dev \
        libz-dev \
        libcurl4-openssl-dev \
        r-base \
        r-base-dev \
        r-cran-rcpp \
        openjdk-11-jdk \
        python3 \
        python3-dev \
        python3-numpy \
        python3-matplotlib \
        python3-pysam \
        python3-htseq \
        samtools \
        littler \
        git \
        curl \
        make \
        unzip \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && apt-get -y update \
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

#install required R packages from CRAN
RUN install.r \
        devtools \
        docopt \
        lattice \
        RCurl \
        openssl \
        Hmisc \
        tidyr \
        readr \
        stringr \
        ggplot2 \
        BiocManager \
    && rm -rf /var/lib/apt/lists/*

#install other required R packages from github
RUN Rscript -e "\
BiocManager::install('annotate'); \
BiocManager::install('GenomicRanges'); \
BiocManager::install('SummarizedExperiment'); \
BiocManager::install('GenomeInfoDb'); \
BiocManager::install('DESeq2'); \
BiocManager::install('BiocParallel')"

#copy the dgexp.R script to /opt/
COPY dgexp.R /opt/

#change working directory
WORKDIR /opt/

#chmod of /opt/
RUN chmod -R 777 /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
