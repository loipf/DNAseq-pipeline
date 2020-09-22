#!/usr/bin/env bash

### maybe set up everything with bioconda ? - dont know if they have newest versions
### export tool-pathes to .bashrc ?


TOOL_DIR="/home/stefan/tools/"

BWA_VERSION="v0.7.17"
CUTADAPT_VERSION="2.10"
SAMTOOLS_VERSION="1.10"
PICARD_VERSION="2.23.4"
FASTQC_VERSION="v0.11.9"
MULTIQC_VERSION="1.9"



cd $TOOL_DIR

### nextflow v20.07.1
curl -fsSL https://get.nextflow.io | bash


### BWA 
git clone https://github.com/lh3/bwa.git --branch $BWA_VERSION
cd bwa; make
cd $TOOL_DIR


### cutadapt
python3 -m pip install --user --upgrade cutadapt==$CUTADAPT_VERSION


### samtools
wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2
tar xvjf samtools-$SAMTOOLS_VERSION.tar.bz2  && rm samtools-$SAMTOOLS_VERSION.tar.bz2 
cd samtools-$SAMTOOLS_VERSION
make
cd $TOOL_DIR


### picard
wget https://github.com/broadinstitute/picard/releases/download/$PICARD_VERSION/picard.jar


### FastQC
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_$FASTQC_VERSION.zip
unzip fastqc_$FASTQC_VERSION.zip && rm fastqc_$FASTQC_VERSION.zip
chmod 755 FastQC/fastqc


### MultiQC
pip install multiqc==$MULTIQC_VERSION










