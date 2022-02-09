#!/usr/bin/env bash

data_dir="../data"

ensembl_release="101"



cd $data_dir


### get reference genome - WGS + WES
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# ### get dbsnp for base recalibration - not needed
# dbsnp_release="154"
# curl https://ftp.ncbi.nih.gov/snp/archive/b$dbsnp_release/VCF/GCF_000001405.38.gz > dbsnp_154_hg38.vcf.gz



### short test genome
#curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz




