

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$projectDir/data"



process CREATE_FOLDER_STRUCTURE {

	input:
		path data_dir

	shell:
	'''
	mkdir -p !{data_dir}
	mkdir -p !{data_dir}/reads_raw
	mkdir -p !{data_dir}/reads_prepro
	'''

}



process DATA_ACQUISITION { 
	storeDir params.data_dir, mode: "copy"

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
	publishDir "$params.data_dir/BWA_index", mode: "copy"

	input:
		path tool_bwa
		path reference_genome

	output:
		path "test_file.txt"
		//path "*.{amb,ann,bwt,pac,sa}", emit: bwa_index_bwamem

	shell:
	'''
	!{tool_bwa} > test_file.txt # index !{reference_genome}
	'''
}
	





process PREPROCESS_READS { 
	publishDir "$data_dir/reads_preprocessed/", mode: "copy", saveAs: { filename -> "${sample_id}/${sample_id}_$filename" }

	input:
		path data_dir
		tuple val(sample_id), path(reads) 
		path tool_cutadapt
		val num_threads
		path adapter_seq


	output:
		tuple path("prepro_1.fastq.gz"), path("prepro_2.fastq.gz"), emit: reads_preprocessed 
		path("cutadapt_output.txt")

	shell:
	'''
	ADAPTER_5=$(cat !{adapter_seq} | sed -n 1p | cut -f 2)  # forward
	ADAPTER_3=$(cat !{adapter_seq} | sed -n 2p | cut -f 2)  # reverse

	!{tool_cutadapt} --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any -b $ADAPTER_5 -B $ADAPTER_3 -o prepro_1.fastq.gz -p prepro_2.fastq.gz !{reads} > cutadapt_output.txt

	'''
}




process MULTIQC_READS { 
	//publishDir "$parent", mode: "copy", saveAs: { filename -> "${sample_id}_$filename" }
	//publishDir "$reads.parent", mode: "copy", saveAs: { filename -> "${sample_id}_$filename" }
	publishDir "$read_dir/$sample_id", mode: "copy", saveAs: { filename -> "${sample_id}_$filename" }

	input:
		path read_dir
		path tool_dir
		val num_threads
		tuple val(sample_id), path(reads) 
		path adapter_seq


	output:
		//tuple path("1_fastqc.zip"), path("2_fastqc.zip")
	path("multiqc_test.txt")


	shell:
	'''
	echo !{read_dir} > multiqc_test.txt
	#!{tool_dir}/FastQC/fastqc -a !{adapter_seq} -t !{num_threads} --noextract !{reads}
	
	'''
}




process MAPPING_BWA { 
	publishDir "$data_dir/reads_mapped/", mode: 'copy'

	input:
		path data_dir
		path tool_dir
		path reference_genome


	output:
		path "Homo_sapiens.GRCh38.dna.alt.fa.gz", emit: reference_genome


	shell:
	"""
	curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz


	"""
}


















process TEST { 
	publishDir "$data_dir/", mode: 'copy'
 
	input:
		path data_dir
		path reference_genome


	output:
		path "test_file.txt"


	shell:
	"""
	head !{reference_genome} > test_file.txt

	"""
}

















