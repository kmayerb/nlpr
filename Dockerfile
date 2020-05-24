FROM ubuntu:18.04

MAINTAINER kmayerbl@fredhutch.org

WORKDIR /

# install basics, g++ compiler and boost components for silix
RUN apt-get update && \
	apt-get -y install wget && \
	apt-get -y install build-essential g++ python-dev autotools-dev libicu-dev build-essential libbz2-dev libboost-all-dev

# download and decompress silix
RUN mkdir /software && \ 
	wget ftp://pbil.univ-lyon1.fr/pub/logiciel/silix/silix-1.2.6.tar.gz -P /software/ && \
	tar -xzvf /software/silix-1.2.6.tar.gz -C /software/

# make install silix
WORKDIR /software/silix-1.2.6/
RUN ./configure && \
  make && \
  make install 

WORKDIR /

# install primer3
RUN apt-get -y install primer3

# install emboss contains fuzznuc
RUN apt-get -y install emboss 

# install recent cd-hit (paper version was incompatable with g++ version: https://github.com/kuleshov/cdhit/issues/12)
RUN apt-get -y install cd-hit

# install ps and text editor helpful for interactive use and AWS batch
RUN apt-get -y install procps && apt-get -y install nano
RUN apt-get -y install git

WORKDIR /