# nlpr

[![Docker Repository on Quay](https://quay.io/repository/kmayerb/nlpr/status "Docker Repository on Quay")](https://quay.io/repository/kmayerb/nlpr)


## Using nLprime

1. If working on OSX or Windows, download, install and run [Docker Desktop](https://www.docker.com/products/docker-desktop). 
2. With Docker Desktop running, from the command line:

```bash
docker pull quay.io/kmayerb/nlpr:0.0.3
docker run -it quay.io/kmayerb/nlpr:0.0.3
git clone https://github.com/kmayerb/nlpr.git
cd /nlpr/nlpr
python nl_control_file_example.py
``` 

## Suggested Primers

After the program runs, you can see example outputs in `nlpr/nlpr/example` folder.


## Inputs

For new runs three input files are requred formatted like the fils in the `nlpr/inputs` folder.

1. Nuleic Acids fasta `*fna`
2. Protein fasta  `*.faa`
3. Tabular results of an all_v_all blastp or blastn `.all-v-all_blastp_output`


## Configuration

Primer design can be configured in `nl_control_file_.py`


## Dockerfile

```bash
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

# install git tools
RUN apt-get -y install git

# install Biopython
RUN apt-get -y install python-pip
RUN pip install Biopython

WORKDIR /
```

### Archive of legacy versions of software used when scripts were originally run in 2014.

* cd-hit-v4.5.4-2011-03-07.tgz [download](https://www.dropbox.com/s/34ybl944fkcefds/cd-hit-v4.5.4-2011-03-07.tgz?dl=1)

* primer3-2.3.4.tar.gz [download]((https://www.dropbox.com/s/z7x7tx1cmvwvl9h/primer3-2.3.4.tar.gz?dl=1)

* ncbi-blast-2.2.18-universal-macosx.tar.gz [download](https://www.dropbox.com/s/y2jeajmxgcho0bt/ncbi-blast-2.2.18%2B-universal-macosx.tar.gz?dl=1)

* silix-1.2.6.tar.gz [download](https://www.dropbox.com/s/rlctg1chfxqr13c/silix-1.2.6.tar.gz?dl=1)
