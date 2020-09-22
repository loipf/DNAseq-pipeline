
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
   } from './modules.nf' 





/*
 * default parameters
 */ 
params.tool_dir   = "/home/stefan/tools"
params.data_dir   = "$projectDir/data"
params.reads   = "$projectDir/data/reads_raw/*/*_{1,2}.fastq.gz"
params.scripts_dir    = "$projectDir/scripts"
params.results_dir    = "$projectDir/results"
params.num_threads = 3



/*
 * other parameters
 */
params.ensembl_release = "101"
params.adapter_seq_file = "$projectDir/data/adapter_seq.txt"




log.info """\
DNA-SEQ PIPELINE
================================
tool_dir : $params.tool_dir
data_dir : $params.data_dir
reads    : $params.reads
results_dir  : $params.results_dir

"""



/* 
 * main pipeline logic
 */
workflow {
    channel_reads = Channel
            .fromFilePairs( params.reads )
            .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }


    DATA_ACQUISITION(params.data_dir, params.ensembl_release)
    PREPROCESS_READS(params.data_dir, params.num_threads, channel_reads, params.adapter_seq_file)


 
    


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








