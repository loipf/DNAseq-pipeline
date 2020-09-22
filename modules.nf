


process DATA_ACQUISITION { 
 
  input:
    path data_dir
    val ensembl_release

   publishDir "$data_dir/", mode: 'copy'

  output:
    path "Homo_sapiens.GRCh38.dna.alt.fa.gz", emit: reference_genome


  shell:
  '''
    curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz

#    curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz

  '''
}




process PREPROCESS_READS { 

  input:
    path data_dir
    val num_threads
    tuple val(sample_id), path(reads) 
    path adapter_seq

  publishDir "$data_dir/reads_preprocessed/", mode: 'copy', saveAs: { filename -> "${sample_id}/${sample_id}_$filename" }

  output:
    tuple path("prepro_1.fastq.gz"), path("prepro_2.fastq.gz"), emit: reads_preprocessed 

  shell:
  '''
    ADAPTER_5=$(cat !{adapter_seq} | sed -n 2p)  # forward
    ADAPTER_3=$(cat !{adapter_seq} | sed -n 4p)  # reverse

    cutadapt --cores=!{num_threads} --max-n 0.1 --discard-trimmed --pair-filter=any -b $ADAPTER_5 -B $ADAPTER_3 -o prepro_1.fastq.gz -p prepro_2.fastq.gz !{reads}
    
  '''
}







process MAPPING_BWA { 
  input:
    path data_dir
    path tool_dir
    path reference_genome

   publishDir "$data_dir/reads_mapped/", mode: 'copy'

  output:
    path "Homo_sapiens.GRCh38.dna.alt.fa.gz", emit: reference_genome


  shell:
  """
    curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz


  """
}


















process TEST { 
 
  input:
    path data_dir
    path reference_genome

  publishDir "$data_dir/", mode: 'copy'

  output:
    file "test_file.txt"


  shell:
  """
    head !{reference_genome} > test_file.txt

  """
}

















