# DNAseq pipeline

a DNAseq mapping pipeline from `.fastq` files to `.bam` files with intermediate quality reports


---
### set up pipeline


before running, you have to set up the attached Docker image:
```sh
sudo docker build -t dna-seq_pipeline https://raw.githubusercontent.com/loipf/DNAseq-pipeline/master/docker/Dockerfile
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker [container hash]` argument.



if you already have an reference genome, you can specify it with the nextflow argument `--reference_genome [genome.fa]` or in the `main.nf` file. otherwise run the following script to download it into `data` folder:
```sh
bash scripts/data_acquisition.sh
```

additional, an adapter sequence file needs to be specified with the nextflow argument `--adapter_seq_file [file.tsv]` or in the `main.nf` file. it must be structured as follows:
```
ADAPTER_5	AANTGG
ADAPTER_3	GATCGG
```


---
### run mapping pipeline

it can be run locally with downloaded github-repo with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAseq-pipeline
```


by default, all output will be saved into the `data` folder




