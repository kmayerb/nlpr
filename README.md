# nlpr


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
RUN apt-get install -y procps && apt-get install -y nano

# try running help commands as a basic install test
RUN /software/silix-1.2.6/src/silix -h && \
	/usr/bin/fuzznuc -h && \
	/usr/bin/cd-hit -h && \
	/usr/bin/primer3_core -h 

WORKDIR /

```


## Additional Notes

### Paths to help docs
```bash
/software/silix-1.2.6/src/silix -h
/usr/bin/fuzznuc -h
/usr/bin/cd-hit -h
/usr/bin/primer3_core -h
```


### Archive of legacy versions of software used 
```
[cd-hit-v4.5.4-2011-03-07.tgz](https://www.dropbox.com/s/34ybl944fkcefds/cd-hit-v4.5.4-2011-03-07.tgz?dl=1)

[primer3-2.3.4.tar.gz](https://www.dropbox.com/s/z7x7tx1cmvwvl9h/primer3-2.3.4.tar.gz?dl=1)

[ncbi-blast-2.2.18-universal-macosx.tar.gz](https://www.dropbox.com/s/y2jeajmxgcho0bt/ncbi-blast-2.2.18%2B-universal-macosx.tar.gz?dl=1)

[silix-1.2.6.tar.gz](https://www.dropbox.com/s/rlctg1chfxqr13c/silix-1.2.6.tar.gz?dl=1)

```
