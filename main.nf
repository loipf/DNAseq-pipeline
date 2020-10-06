
/* 
 * DNA-SEQ PIPELINE 
 * for paired-end reads
 */


/* 
 * import modules 
 */
nextflow.enable.dsl=2

include { 
	CREATE_FOLDER_STRUCTURE;
	DATA_ACQUISITION;
	CREATE_BWA_INDEX;
	TEST;
	PREPROCESS_READS;
	MULTIQC_READS
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev = true

params.reads		= "$projectDir/test_reads_dir/*/*_{1,2}.fastq.gz"
params.data_dir		= "$projectDir/data"
params.scripts_dir	= "$projectDir/scripts"


/*
 * tool paths
 */ 
params.tool_fastqc		= "/home/stefan/FastQC/fastqc"
params.tool_cutadapt	= "/home/stefan/.local/bin/cutadapt"
params.tool_multiqc		= "/home/stefan/miniconda3/bin/multiqc"
params.tool_samtools	= "/home/stefan/tools/samtools-1.10/samtools"
params.tool_bwa 		= "/home/stefan/tools/bwa/bwa"
params.tool_deeptools	= "/home/stefan/miniconda3/bin/deeptools"


/*
 * other parameters
 */
params.num_threads		= 3
params.ensembl_release	= "101"
params.adapter_seq_file	= "$projectDir/data/adapter_seq.tsv"
params.reference_genome	= "$projectDir/data/Homo_sapiens.GRCh38.dna.alt.fa.gz"




log.info """\
DNA-SEQ PIPELINE
================================
reads		: $params.reads
data_dir	: $params.data_dir

"""



/* 
 * main pipeline logic
 */
workflow {
	channel_reads = Channel
			.fromFilePairs( params.reads )
			.ifEmpty { error "cannot find any reads matching: ${params.reads}" }
			.take( params.dev ? 2 : -1 )  // only consider 2 files for debugging


	CREATE_FOLDER_STRUCTURE(params.data_dir)

	//DATA_ACQUISITION(params.data_dir, params.ensembl_release)  # STOREDIR DOES NOT WORK

	CREATE_BWA_INDEX(params.tool_bwa, params.reference_genome)

	//PREPROCESS_READS(params.data_dir+"/reads_test1/", channel_reads, params.tool_cutadapt, params.num_threads, params.adapter_seq_file)

	//MULTIQC_READS(params.tool_dir, params.num_threads, PREPROCESS_READS.out.reads_preprocessed, params.adapter_seq_file)

	//MULTIQC_READS( (params.data_dir+"/reads_raw/"), params.tool_dir, params.num_threads, channel_reads, params.adapter_seq_file)



}


  //  TEST(params.data_dir, DATA_ACQUISITION.out.reference_genome)

/*
* params.dev = false
* params.number_of_inputs = 2
* Channel
*	.from(1..300)
*	.take( params.dev ? params.number_of_inputs : -1 )
*	.println() 
*/








