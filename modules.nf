

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"



// maybe outsource this script to simple bash 
process DATA_ACQUISITION { 
	storeDir params.data_dir, mode: "copy"  // DOES NOT WORK

	input:
		path data_dir
		val ensembl_release

	output:
		path "Homo_sapiens.GRCh38.dna.alt.fa.gz", emit: reference_genome


	shell:
	'''
#	curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

	curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

	'''
}



process PREPROCESS_READS { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", pattern:"*cutadapt_output.txt", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_3_seq_file
		path adapter_5_seq_file

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''

	### adapter 3' input as string or file
	if [[ !{adapter_3_seq_file} == *".fasta"* ]]
	then
  		ADAPTER_3=file:!{adapter_3_seq_file}
	else
		ADAPTER_3="!{adapter_3_seq_file}"
	fi

	### adapter 5' input as string or file
	if [[ !{adapter_5_seq_file} == *".fasta"* ]]
	then
  		ADAPTER_5=file:!{adapter_5_seq_file}
	else
		ADAPTER_5="!{adapter_5_seq_file}"
	fi
	
	### theoretically sorting not needed but maybe speeds up things
	reads_sorted=$(echo !{reads} | xargs -n1 | sort | xargs)

	### combine multiple seq files in the same sample directory with same direction together
	reads_sorted_1=$(find $reads_sorted -name "*_1.fq.gz" -o -name "*_1.fastq.gz")
	reads_sorted_2=$(find $reads_sorted -name "*_2.fq.gz" -o -name "*_2.fastq.gz")
	cat $reads_sorted_1 > raw_reads_connected_1.fastq.gz
	cat $reads_sorted_2 > raw_reads_connected_2.fastq.gz

	cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -a $ADAPTER_3 -A $ADAPTER_5 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz raw_reads_connected_1.fastq.gz raw_reads_connected_2.fastq.gz > !{sample_id}_cutadapt_output.txt
	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }
	stageInMode = 'copy'   // avoids permission denied error
	cache false

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_3_seq_file
		path adapter_5_seq_file

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	cat !{adapter_3_seq_file} !{adapter_5_seq_file} > adapter_both.fasta
	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' adapter_both.fasta > adapter_both.txt

	fastqc -a adapter_both.txt -t !{num_threads} --noextract !{reads}
	'''
}



process FASTQC_READS_PREPRO { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path adapter_3_seq_file
		path adapter_5_seq_file

	output:
		path "*.zip", emit: reports
		path "*.html"

	shell:
	'''
	cat !{adapter_3_seq_file} !{adapter_5_seq_file} > adapter_both.fasta
	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' adapter_both.fasta > adapter_both.txt

	fastqc -a adapter_both.txt -t !{num_threads} --noextract !{reads_prepro}
	'''
}



process CREATE_BWA_INDEX { 
	publishDir "$params.data_dir/bwa_index", mode: "copy"

	input:
		path reference_genome

	output:
		path "*.{amb,ann,bwt,pac,sa}", emit: bwa_index

	shell:
	'''
	bwa index !{reference_genome} 
	'''
}



process MAPPING_BWA { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_mapped", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(reads_prepro) 
		val num_threads
		path reference_genome
		path bwa_index  // to ensure index is created


	output:
		path "${sample_id}.bam", emit: reads_mapped
		path "${sample_id}.bam.bai", emit: reads_mapped_index
		path "*", emit: all


	shell:
	'''
	bwa mem -Y -R "@RG\\tID:!{sample_id}\\tSM:!{sample_id}" -t !{num_threads} -K 100000000 !{reference_genome} !{reads_prepro} \
	| samtools view -@ !{num_threads} -h -b - \
	| samtools sort -n -@ !{num_threads} - \
	| samtools fixmate -m -@ !{num_threads} - - \
	| samtools sort -@ !{num_threads} - \
	| samtools markdup -@ !{num_threads} -f !{sample_id}_markdup_stats.txt - !{sample_id}.bam

	samtools index -b -@ !{num_threads} !{sample_id}.bam
	samtools stats -@ !{num_threads} !{sample_id}.bam > !{sample_id}_stats.txt

	'''
}




process DEEPTOOLS_ANALYSIS { 
	publishDir "$params.data_dir/reads_mapped/_deepTools", mode: 'copy'

	input:
		path reads_mapped
		path reads_mapped_index
		val num_threads

	output:
		path "*", emit:all


	shell:
	'''
	multiBamSummary bins -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} -o multiBamSummary.npz

	plotCorrelation --corData multiBamSummary.npz --corMethod spearman --whatToPlot heatmap --outFileCorMatrix plotCorrelation_matrix.tsv

	plotPCA --corData multiBamSummary.npz --outFileNameData plotPCA_matrix.tsv

	plotCoverage -p !{num_threads} --ignoreDuplicates --smartLabels --bamfiles !{reads_mapped} --outRawCounts plotCoverage_rawCounts_woDuplicates.tsv > plotCoverage_output.tsv

	bamPEFragmentSize -p !{num_threads} --bamfiles !{reads_mapped} --table bamPEFragment_table.tsv --outRawFragmentLengths bamPEFragment_rawLength.tsv

	estimateReadFiltering -p !{num_threads} --smartLabels --bamfiles !{reads_mapped} > estimateReadFiltering_output.tsv

	'''
}






process MULTIQC_RAW { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_raw .
	'''
}


process MULTIQC_PREPRO { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_prepro .
	'''
}


process MULTIQC_MAPPED { 
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o reads_mapped .
	'''
}


















