#!/bin/bash

############################################
### specifiy directories, tools and properties

reads_dir="/home/stefan/Documents/umcg/dna-seq_pipeline/test_reads_dir"
out_dir="/home/stefan/Documents/umcg/bash_dna-seq_pipeline"

tool_fastqc="/home/stefan/tools/FastQC/fastqc"
tool_cutadapt="/home/stefan/.local/bin/cutadapt"
tool_multiqc="/home/stefan/miniconda3/bin/multiqc"
tool_samtools="/home/stefan/tools/samtools-1.10/samtools"
tool_bwa="/home/stefan/tools/bwa/bwa"
tool_deeptools="/home/stefan/miniconda3/bin/deeptools"


ensembl_release="101"
num_threads="3"
adapter_seq_file="/home/stefan/Documents/umcg/bash_dna-seq_pipeline/data/adapter_seq.tsv"



############################################

mkdir -p $out_dir/data
mkdir -p $out_dir/data/reads_raw
mkdir -p $out_dir/data/reads_prepro


### download reference genome
curl ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > $out_dir/data/Homo_sapiens.GRCh38.dna.alt.fa.gz


### create BWA index and move in separate folder
$tool_bwa index $out_dir/data/Homo_sapiens.GRCh38.dna.alt.fa.gz
cd $out_dir/data
mkdir -p BWA_index
mv -t BWA_index *.amb *.ann *.bwt *.pac *.sa


### trim reads
adapter_5=$(cat $adapter_seq_file | sed -n 1p | cut -f 2)  # forward
adapter_3=$(cat $adapter_seq_file | sed -n 2p | cut -f 2)  # reverse

############################################
### iterate over samples and map

for dir_path in $reads_dir/*
do
#    dir_path="/home/stefan/Documents/umcg/dna-seq_pipeline/test_reads_dir/SRR037207"
    sample_id=$(basename $dir_path)
    read_files=$(find $dir_path -name "*.fastq.gz")

    ### check raw reads
    mkdir -p $out_dir/data/reads_raw/$sample_id
    cd $out_dir/data/reads_raw/$sample_id

    $tool_fastqc -a $adapter_seq_file -t $num_threads --outdir . --noextract $read_files 


    ### preprocess reads
    mkdir -p $out_dir/data/reads_prepro/$sample_id
    cd $out_dir/data/reads_prepro/$sample_id
    
    ### TODO change minimum length - or remove low quality reads
    $tool_cutadapt --cores=$num_threads --max-n 0.1 --minimum-length 10 --discard-trimmed --pair-filter=any -b $adapter_5 -B $adapter_3 -o $sample_id\_prepro_1.fastq.gz -p $sample_id\_prepro_2.fastq.gz $read_files > cutadapt_output.txt

    read_files_prepro=$(find $PWD -name "*.fastq.gz")
    $tool_fastqc -a $adapter_seq_file -t $num_threads --outdir . --noextract $read_files_prepro


    ### mapping process
    ### TODO remove duplicates markdup -r?
    ### TODO disable proper pair check fixmate -p ?
	### TODO check multiple threads add up ? 
    mkdir -p $out_dir/data/reads_mapped/$sample_id
    cd $out_dir/data/reads_mapped/$sample_id
    $tool_bwa mem -Y -t $num_threads -K 100000000 $out_dir/data/BWA_index/Homo_sapiens.GRCh38.dna.alt.fa.gz $read_files_prepro \
	| $tool_samtools view -@ $num_threads -h -b - \
    | $tool_samtools sort -n -@ $num_threads - \
	| $tool_samtools fixmate -m -@ $num_threads - - \
	| $tool_samtools sort -@ $num_threads - \
	| $tool_samtools markdup -@ $num_threads -f $sample_id\_markdup_stats.txt - $sample_id.bam

	$tool_samtools index -@ $num_threads $sample_id.bam
	$tool_samtools stats -@ $num_threads $sample_id.bam > $sample_id\_stats.txt

### TODO discuss, keep only main chromosomes ? 
# https://www.biostars.org/p/361067/
#samtools idxstats in.bam | cut -f 1 | grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
#xargs samtools view -f 1 -F 1284 -q 20 -o out.bam in.bam


done



############################################
### create quality reports

mkdir -p $out_dir/data/quality_reports

read_files_mapped=$(find $out_dir/data/reads_mapped -name "*.bam")


############################################
### deeptools analysis
mkdir -p $out_dir/data/reads_mapped/_deepTools
multiBamSummary bins -p $num_threads --smartLabels --bamfiles $read_files_mapped -o $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz
plotCorrelation --corData $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --outFileCorMatrix $out_dir/data/reads_mapped/_deepTools/plotCorrelation_matrix.tsv
plotPCA --corData $out_dir/data/reads_mapped/_deepTools/multiBamSummary.npz --outFileNameData $out_dir/data/reads_mapped/_deepTools/plotPCA_matrix.tsv

plotCoverage -p $num_threads --ignoreDuplicates --bamfiles $read_files_mapped --outRawCounts $out_dir/data/reads_mapped/_deepTools/plotCoverage_rawCounts_woDuplicates.tsv > $out_dir/data/reads_mapped/_deepTools/plotCoverage_output.tsv

bamPEFragmentSize -p $num_threads --bamfiles $read_files_mapped --table $out_dir/data/reads_mapped/_deepTools/bamPEFragment_table.tsv --outRawFragmentLengths $out_dir/data/reads_mapped/_deepTools/bamPEFragment_rawLength.tsv

estimateReadFiltering -p $num_threads --smartLabels --bamfiles $read_files_mapped > $out_dir/data/reads_mapped/_deepTools/estimateReadFiltering_output.tsv

### need BED files
plotEnrichment -p $num_threads --smartLabels --bamfiles $read_files_mapped --outRawCounts $out_dir/data/reads_mapped/_deepTools/plotEnrichment_rawCounts.tsv

multiqc -f -o $out_dir/data/quality_reports/reads_raw $out_dir/data/reads_raw
multiqc -f -o $out_dir/data/quality_reports/reads_prepro $out_dir/data/reads_prepro
multiqc -f -o $out_dir/data/quality_reports/reads_mapped $out_dir/data/reads_mapped


















