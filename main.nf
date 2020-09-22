
/* 
 * DNA-SEQ PIPELINE
 */


/* 
 * import modules 
 */
nextflow.enable.dsl=2

include { 
  DATA_ACQUISITION;
   } from './modules.nf' 





/*
 * default parameters
 */ 
params.tool_dir   = "/home/stefan/tools"
params.data_dir   = "$projectDir/data"
params.scripts_dir    = "$projectDir/scripts"
params.results_dir    = "$projectDir/results"
params.num_threads = 3



/*
 * other parameters
 */
params.ensembl_release = "101"





log.info """\
DNA-SEQ PIPELINE
================================
tool_dir : $params.tool_dir
data_dir : $params.data_dir
results_dir  : $params.results_dir

"""



/* 
 * main pipeline logic
 */
workflow {
      //reads_ch = Channel.fromFilePairs(params.reads)

      DATA_ACQUISITION(params.data_dir, params.ensembl_release)





}




/*
* params.dev = false
* params.number_of_inputs = 2
* Channel
*    .from(1..300)
*    .take( params.dev ? params.number_of_inputs : -1 )
*    .println() 
*/








