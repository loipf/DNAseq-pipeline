
/* 
 * DNA-SEQ PIPELINE
 * inspired by https://github.com/CRG-CNAG/CalliNGS-NF/blob/master/main.nf
 */



/*
 * Define the default parameters
 */ 
params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.denylist   = "$baseDir/data/denylist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"

log.info """\
DNA-SEQ PIPELINE
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
denylist : $params.denylist
results  : $params.results
"""

/* 
 * Import modules 
 */
include { 
  PREPARE_GENOME_SAMTOOLS;
  PREPARE_GENOME_PICARD; 
  PREPARE_STAR_GENOME_INDEX;
  PREPARE_VCF_FILE;
  RNASEQ_MAPPING_STAR;
  RNASEQ_GATK_SPLITNCIGAR; 
  RNASEQ_GATK_RECALIBRATE;
  RNASEQ_CALL_VARIANTS;
  POST_PROCESS_VCF;
  PREPARE_VCF_FOR_ASE;
  ASE_KNOWNSNPS;
  group_per_sample } from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {
      reads_ch = Channel.fromFilePairs(params.reads)

      // PART 1: Data preparation
      PREPARE_GENOME_SAMTOOLS(params.genome)
      PREPARE_GENOME_PICARD(params.genome)
      PREPARE_STAR_GENOME_INDEX(params.genome)
      PREPARE_VCF_FILE(params.variants, params.denylist)

      // PART 2: STAR RNA-Seq Mapping
      RNASEQ_MAPPING_STAR( 
            params.genome, 
            PREPARE_STAR_GENOME_INDEX.out, 
            reads_ch )

      // PART 3: GATK Prepare Mapped Reads
      RNASEQ_GATK_SPLITNCIGAR(
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_MAPPING_STAR.out )

      // PART 4: GATK Base Quality Score Recalibration Workflow
      RNASEQ_GATK_RECALIBRATE(
                  params.genome, 
                  PREPARE_GENOME_SAMTOOLS.out, 
                  PREPARE_GENOME_PICARD.out, 
                  RNASEQ_GATK_SPLITNCIGAR.out, 
                  PREPARE_VCF_FILE.out)

      // PART 5: GATK Variant Calling
      RNASEQ_CALL_VARIANTS( 
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            RNASEQ_GATK_RECALIBRATE.out.groupTuple() )

      // PART 6: Post-process variants file and prepare for 
      // Allele-Specific Expression and RNA Editing Analysis
      POST_PROCESS_VCF( 
            RNASEQ_CALL_VARIANTS.out, 
            PREPARE_VCF_FILE.out )

      PREPARE_VCF_FOR_ASE( POST_PROCESS_VCF.out )

      ASE_KNOWNSNPS(
            params.genome, 
            PREPARE_GENOME_SAMTOOLS.out, 
            PREPARE_GENOME_PICARD.out, 
            group_per_sample(
                  RNASEQ_GATK_RECALIBRATE.out, 
                  PREPARE_VCF_FOR_ASE.out[0]) )
}
