


process DATA_ACQUISITION { 
 
  input:
    path data_dir
    val ensembl_release

  publishDir "$data_dir/", mode: 'copy'

  output:
    path "Homo_sapiens.GRCh38.dna.alt.fa.gz", emit: reference_genome


  shell:
  """
    curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz

#    curl ftp://ftp.ensembl.org/pub/release-!{ensembl_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz > Homo_sapiens.GRCh38.dna.alt.fa.gz


  """
}






