

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$projectDir/data"



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
	curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz

#	  curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz

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
	


process PREPROCESS_READS { 
	tag "cutadapt on $sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_seq

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "cutadapt_output.txt"

	shell:
	'''
	ADAPTER_5=$(cat !{adapter_seq} | sed -n 1p | cut -f 2)  # forward
	ADAPTER_3=$(cat !{adapter_seq} | sed -n 2p | cut -f 2)  # reverse

	cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any -b $ADAPTER_5 -B $ADAPTER_3 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz !{reads} > cutadapt_output.txt

	'''
}



// not possible to run dynamically fastqc with same name
process FASTQC_READS_RAW { 
	tag "fastqc_raw on $sample_id"
	publishDir "$params.data_dir/reads_raw", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_seq

	output:
		path "*.zip"
		path "*.html"

	shell:
	'''
	fastqc -a !{adapter_seq} -t !{num_threads} --noextract !{reads}
	'''
}



process FASTQC_READS_PREPRO { 
	tag "fastqc_prepro on $sample_id"
	publishDir "$params.data_dir/reads_prepro", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path adapter_seq

	output:
		path "*.zip"
		path "*.html"

	shell:
	'''
	fastqc -a !{adapter_seq} -t !{num_threads} --noextract !{reads}
	'''
}




process MAPPING_BWA { 
	tag "bwa mapping on $sample_id"
	publishDir "$params.data_dir/reads_mapped", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(reads) 
		val num_threads
		path reference_genome
		path bwa_index  // to ensure index is created


	output:
		path "${sample_id}.bam", emit: reads_mapped
		path "${sample_id}_stats.txt"
		path "*"
		//path "${sample_id}_markup_stats.txt"  // problems saving
		//path "${sample_id}.bam.ai"


	shell:
	'''
	bwa mem -Y -t !{num_threads} -K 100000000 !{reference_genome} !{reads} \
	| samtools view -@ !{num_threads} -h -b - \
    | samtools sort -n -@ !{num_threads} - \
	| samtools fixmate -m -@ !{num_threads} - - \
	| samtools sort -@ !{num_threads} - \
	| samtools markdup -@ !{num_threads} -f !{sample_id}_markdup_stats.txt - !{sample_id}.bam

	samtools index -b -@ !{num_threads} !{sample_id}.bam
	samtools stats -@ !{num_threads} !{sample_id}.bam > !{sample_id}_stats.txt

	'''
}
















//.map{ file -> tuple(file[0], file)}

process TEST { 
	publishDir "$params.data_dir/test", mode: 'copy'
 
	input:
		tuple val(sample_id), path(reads) 

	output:
		tuple val(sample_id), path("test_1.fastq.txt"), path("test_2.fastq.txt"), emit: reads_prepro


	shell:
	"""
	echo !{sample_id} > test_1.fastq.txt
	echo !{sample_id} > test_2.fastq.txt
	"""
}


process TEST2 { 
	publishDir "$params.data_dir/test2", mode: 'copy'
 
	input:
		tuple val(sample_id), path(reads) 

	output:
		path "${sample_id}" 
		//path("test2_1.fastq.txt")
		//path("test2_2.fastq.txt")

	shell:
	"""
	echo !{sample_id} > test2_1.fastq.txt
	echo !{reads} > test2_2.fastq.txt
	"""
}











