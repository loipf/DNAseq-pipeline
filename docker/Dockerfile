FROM conda/miniconda3@sha256:2141bfc4b60cc5eb332402207c80d884daa72fcefa48b50f3ccadd934d1f3d03

MAINTAINER loipf

#RUN apt-get update 
#RUN conda update --all

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
RUN conda install \
		bwa==0.7.17 \
		fastqc==0.11.9-0 \
		deeptools==3.5.0-0 \
		multiqc==1.9 \
		samtools==1.11.0 \
		cutadapt==2.10





