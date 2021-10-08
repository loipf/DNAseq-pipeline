#!/usr/bin/env bash

data_dir="../data"

ensembl_release="101"




cd $data_dir


### get reference genome - WGS
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

### get reference genome - WES
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz > Homo_sapiens.GRCh38.cdna.all.fa.gz



### short test genome
#curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

