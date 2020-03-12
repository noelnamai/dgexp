#RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome
#Bo Li and Colin N Dewey
#University of Wisconsin-Madison
#Department of Computer Sciences 
#Madison, WI

#specifying the base image
FROM ubuntu:19.04

#add metadata to image e.g. name and email of the person who maintains the file
LABEL Tool="RSEM"
LABEL Version="1.3.1"
LABEL Maintainer="Noel Namai"
LABEL Email="noelnamai@yahoo.com"

#set system arguments and enviromental variables
ARG DEBIAN_FRONTEND=noninteractive

#http://bugs.python.org/issue19846: setting "LANG=C" on a Linux system *fundamentally breaks Python 3*
ENV LANG C.UTF-8

#updating ubuntu and installing other necessary software
RUN apt-get update --fix-missing \
	&& apt-get install -y \
		g++ \
		wget \
		bzip2 \
		unzip \
		ca-certificates \
		build-essential \
		checkinstall \
		zlib1g-dev \
		libtbb-dev \
		libncurses5-dev \
		python \
	&& apt-get clean

#change working directory
WORKDIR /usr/local/bin/

#install star
RUN wget https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz \
	&& tar -xzf 2.7.1a.tar.gz \
	&& rm 2.7.1a.tar.gz \
	&& cd STAR-2.7.1a/source \
	&& make STAR

#add STAR dir to path variable
ENV PATH="/usr/local/bin/STAR-2.7.1a/bin/Linux_x86_64:${PATH}"

#install bowtie
RUN wget https://github.com/BenLangmead/bowtie/archive/v1.2.3.tar.gz \
	&& tar -xzf v1.2.3.tar.gz \
	&& rm v1.2.3.tar.gz \
	&& cd bowtie-1.2.3 \
	&& make

#add bowtie dir to path variable
ENV PATH="/usr/local/bin/bowtie-1.2.3:${PATH}"

#install bowtie2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-source.zip \
	&& unzip bowtie2-2.3.5.1-source.zip \
	&& rm bowtie2-2.3.5.1-source.zip \
	&& cd bowtie2-2.3.5.1 \
	&& make

#add bowtie2 dir to path variable
ENV PATH="/usr/local/bin/bowtie2-2.3.5.1:${PATH}"

#install RSEM
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz \
	&& tar xzvf v1.3.1.tar.gz \
	&& rm v1.3.1.tar.gz \
	&& cd RSEM-1.3.1 \
	&& make

#add RSEM dir to path variable
ENV PATH="/usr/local/bin/RSEM-1.3.1:${PATH}"

#chmod of /home/user and change working directory
RUN chmod -R 777 /opt/

#change working directory
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]
