
/* 
 * DNA-SEQ PIPELINE 
 * for paired-end reads
 */


/* 
 * import modules 
 */
nextflow.enable.dsl=2

include { 
    DATA_ACQUISITION;
    TEST;
    PREPROCESS_READS;
	MULTIQC_READS
} from './modules.nf' 





/*
 * default parameters
 */ 
params.reads       = "$projectDir/test_reads_dir/*/*_{1,2}.fastq.gz"
params.data_dir    = "$projectDir/data"
params.scripts_dir = "$projectDir/scripts"
params.results_dir = "$projectDir/results"


/*
 * tool paths
 */ 
params.tool_fastqc   = "/home/stefan/FastQC/fastqc"
params.tool_cutadapt = "/home/stefan/.local/bin/cutadapt"
params.tool_picard   = "/home/stefan/tools/picard.jar"
params.tool_multiqc  = "/home/stefan/miniconda3/bin/multiqc"
params.tool_samtools = "/home/stefan/tools/samtools-1.10/samtools"
params.tool_bwa      = "/home/stefan/tools/bwa"


/*
 * other parameters
 */
params.num_threads      = 3
params.ensembl_release  = "101"
params.adapter_seq_file = "$projectDir/data/adapter_seq.tsv"




log.info """\
DNA-SEQ PIPELINE
================================
data_dir    : $params.data_dir
reads       : $params.reads
results_dir : $params.results_dir

"""



/* 
 * main pipeline logic
 */
workflow {
    channel_reads = Channel
            .fromFilePairs( params.reads )
            .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }


	DATA_ACQUISITION(params.data_dir, params.ensembl_release)
	PREPROCESS_READS(params.data_dir+"/reads_test1/", channel_reads, params.tool_cutadapt, params.num_threads, params.adapter_seq_file)
	//MULTIQC_READS(params.tool_dir, params.num_threads, PREPROCESS_READS.out.reads_preprocessed, params.adapter_seq_file)

	//MULTIQC_READS( (params.data_dir+"/reads_raw/"), params.tool_dir, params.num_threads, channel_reads, params.adapter_seq_file)



}


  //  TEST(params.data_dir, DATA_ACQUISITION.out.reference_genome)

/*
* params.dev = false
* params.number_of_inputs = 2
* Channel
*    .from(1..300)
*    .take( params.dev ? params.number_of_inputs : -1 )
*    .println() 
*/








