
/* 
 * DNAseq-pipeline
 * for paired-end reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	DATA_ACQUISITION;
	CREATE_BWA_INDEX;
	PREPROCESS_READS;
	FASTQC_READS_RAW;
	FASTQC_READS_PREPRO;
	MAPPING_BWA;
	DEEPTOOLS_ANALYSIS;
	MULTIQC_RAW;
	MULTIQC_PREPRO;
	MULTIQC_MAPPED;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = -1

params.project_dir	= "$projectDir"
params.reads_dir	= "$params.project_dir/data/reads_raw"

params.reads		= "$params.reads_dir/*/*_{1,2}.{fastq,fq}.gz"
params.data_dir		= "$params.project_dir/data"
params.scripts_dir	= "$params.project_dir/scripts"


/*
 * other parameters
 */

params.num_threads		= 3
params.adapter_3_seq_file	= "$params.project_dir/data/adapter_3_seq.fa"
params.adapter_5_seq_file	= "$params.project_dir/data/adapter_5_seq.fa"
params.reference_genome	= "$params.project_dir/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
params.seq_type = "wgs" // wgs or wes




log.info """\
DNAseq-pipeline
===================================================
reads			: $params.reads
data_dir		: $params.data_dir
reference_genome	: $params.reference_genome
adapter_3_seq_file	: $params.adapter_3_seq_file
adapter_5_seq_file	: $params.adapter_5_seq_file

===================================================

"""



/* 
 * main pipeline logic
 */
workflow {

	// extension from .fromFilePairs( params.reads ) to deal with multiple sequencing files which get combined in module
	channel_reads = Channel
			.fromPath(params.reads)
			.map{ files -> tuple(files.getParent().getName(), files) }
			.groupTuple()
			.ifEmpty { error "cannot find any reads matching: ${params.reads}" }
			.take( params.dev_samples )  // only consider a few files for debugging

	channel_reads.view()

	// // DATA_ACQUISITION(params.data_dir, params.ensembl_release)  # STOREDIR DOES NOT WORK
	//CREATE_BWA_INDEX(params.reference_genome)
	//PREPROCESS_READS(channel_reads, params.num_threads, params.adapter_3_seq_file, params.adapter_5_seq_file)
	//channel_reads_prepro = PREPROCESS_READS.out.reads_prepro.map{ it -> tuple(it[0], tuple(it[1], it[2])) }

	//FASTQC_READS_RAW(channel_reads, params.num_threads, params.adapter_seq_file)
	//FASTQC_READS_PREPRO(channel_reads_prepro, params.num_threads, params.adapter_seq_file)

	//MAPPING_BWA(channel_reads_prepro, params.num_threads, params.reference_genome, CREATE_BWA_INDEX.out.bwa_index.collect())

	//DEEPTOOLS_ANALYSIS(MAPPING_BWA.out.reads_mapped.collect(), MAPPING_BWA.out.reads_mapped_index.collect(), params.num_threads)

	//MULTIQC_RAW(FASTQC_READS_RAW.out.reports.collect() )
	//MULTIQC_PREPRO(FASTQC_READS_PREPRO.out.reports.concat(PREPROCESS_READS.out.cutadapt).collect() )
	//MULTIQC_MAPPED(MAPPING_BWA.out.all.concat(DEEPTOOLS_ANALYSIS.out.all).collect())
	

}




workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 

















