# DNAseq pipeline

a DNAseq mapping pipeline from `.fastq` files to `.bam` files with intermediate quality reports


---
### set up pipeline


before running, you have to set up the attached Docker image:
```sh
docker build -t dnaseq-pipeline https://raw.githubusercontent.com/loipf/DNAseq-pipeline/master/docker/Dockerfile
```

now either replace the Docker container hash (last output line from previous build command) in `nextflow.config` or run nextflow with the `-with-docker dnaseq-pipeline` argument.



if you already have an reference genome, you can specify it with the nextflow argument `--reference_genome [genome.fa]` or in the `main.nf` file. otherwise run the following script to download it into `data` folder:
```sh
bash scripts/data_acquisition.sh
```

additional, an 3' and 5' adapter sequence (file) needs to be specified with the nextflow arguments `--adapter_3_seq_file [sequence|file.fasta]` and `--adapter_5_seq_file [sequence|file.fasta]` or in the `main.nf` file. if a file is provided, it must be structured like the following example:
```
> adapter_3_batch_01
AANTGG
> adapter_3_batch_02
GATCGG
```


---
### run mapping pipeline

it can be run locally with downloaded github-repo and edited `nextflow.config` file with:
```sh
nextflow run main.nf
```

or

```sh
nextflow run loipf/DNAseq-pipeline -r main --project_dir /path/to/folder --reads_dir /path/to/samples --num_threads 10 --adapter_3_seq_file adapter_3.fasta --adapter_5_seq_file adapter_5.fasta --reference_genome genome.fasta -with-docker dnaseq-pipeline
```
for this execution to work properly, you have to be in the current project directory.


optional extendable with:
```sh
-resume
-with-report report_DNAseq-pipeline
-with-timeline timeline_DNAseq-pipeline
-w work_dir
```


by default, all output will be saved into the `data` folder of the current directory




