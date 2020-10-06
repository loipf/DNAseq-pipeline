FROM conda/miniconda3:latest

MAINTAINER loipf

RUN apt-get update && apt-get install -y curl 

RUN conda update --all
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

### nextflow
RUN curl -fsSL https://get.nextflow.io | bash


# sudo docker build .
# sudo docker run --name test ad0077dfbe3e
# sudo docker stop test
# sudo docker rm test
# sudo docker ps



