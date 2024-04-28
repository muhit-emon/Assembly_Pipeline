#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.R1 = "Read1_file"
params.R2 = "Read2_file"

params.out_fname = "output"

//params.c = "contig_file"
//params.g = "gene_file"

params.ARG_DB = "$projectDir/DB/DeepARG-DB"
params.mobileOG_DB = "$projectDir/DB/mobileOG-DB"

process QC{

  //publishDir "$projectDir", mode: "copy"

  input:
  path R1
  path R2

  output:
  path "${params.out_fname}_qc_R1.fastq.gz", emit: qc_R1
  path "${params.out_fname}_qc_R2.fastq.gz", emit: qc_R2

  """
  $projectDir/fastp -i $R1 -I $R2 -o ${params.out_fname}_qc_R1.fastq.gz -O ${params.out_fname}_qc_R2.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --low_complexity_filter --average_qual 10 --thread 8
  """

}

process assembly{
	
	//publishDir "$projectDir", mode: "copy"
	
  input:
  path qc_R1
  path qc_R2
  
  val output_fname
  
  output:
  path "${output_fname}_assembly/${output_fname}.contigs.fa", emit: contigs

  """
  megahit --presets meta-sensitive -1 $qc_R1 -2 $qc_R2 -t 8 -m 0.5 -f -o ${output_fname}_assembly --out-prefix ${output_fname}
  """

}

process gene_prediction{
	
	//publishDir "$projectDir", mode: "copy"
	
  input:
  path contigs
  
  output:
  path "${params.out_fname}_predicted_genes.fna", emit: predicted_genes
  path "${params.out_fname}_predicted_prots.faa", emit: predicted_prots

  """
  $projectDir/prodigal -i $contigs -d ${params.out_fname}_predicted_genes.fna -a ${params.out_fname}_predicted_prots.faa -p meta
  """

}

process align_full_length_genes_with_ARG_DB{
	
	//publishDir "$projectDir", mode: "copy"
	
  input:
  path predicted_prots
  
  output:
  path "${params.out_fname}_aln_ARG.tsv", emit: aln_ARG_tsv

  """
  $projectDir/diamond blastp -q $predicted_prots -d ${params.ARG_DB} -o ${params.out_fname}_aln_ARG.tsv -k 1 -p 8 --quiet --id 60 --query-cover 70 --subject-cover 60 --evalue 1e-10 -f 6 qtitle stitle pident bitscore evalue
  """

}

process align_full_length_genes_with_mobileOG_DB{
	
	//publishDir "$projectDir", mode: "copy"
	
  input:
  path predicted_prots
  
  output:
  path "${params.out_fname}_aln_mobileOG.tsv", emit: aln_mobileOG_tsv

  """
  $projectDir/diamond blastp -q $predicted_prots -d ${params.mobileOG_DB} -o ${params.out_fname}_aln_mobileOG.tsv -k 1 -p 8 --quiet --id 60 --query-cover 70 --subject-cover 60 --evalue 1e-10 -f 6 qtitle stitle pident bitscore evalue
  """

}

process metacompare{

  publishDir "$projectDir", mode: "copy"

  input:
  path ctg
  path gene
  val output

  output:
  path "${output}_resistome_risk.txt"

  """
  python3 $projectDir/metacmp_for_nf.py -c $ctg -g $gene -o ${output}
  """

}

process output_ARGs_and_mobility{

  publishDir "$projectDir", mode: "copy"

  input:
  path predicted_prots
  path arg_aF
  path mobileOG_aF

  output:
  path "${params.out_fname}_ARGs.faa"
  path "${params.out_fname}_ARGs_and_mobility.tsv"

  """
  python3 $projectDir/extract_ARG_like_seqs_and_their_mobility.py -pG $predicted_prots -ARG_aF $arg_aF -mobileOG_aF $mobileOG_aF -out ${params.out_fname}
  """

}


workflow {

  /*
  contig_ch = Channel.from(params.c)
  
  gene_file_ch = gene_prediction(contig_ch)

  arg_ch = align_full_length_genes_with_ARG_DB(gene_file_ch.predicted_prots)

  mobileOG_ch = align_full_length_genes_with_mobileOG_DB(gene_file_ch.predicted_prots)

  metacompare(contig_ch, gene_file_ch.predicted_genes, params.out_fname)

  output_ARGs_and_mobility(gene_file_ch.predicted_prots, arg_ch.aln_ARG_tsv, mobileOG_ch.aln_mobileOG_tsv)*/

  
  fw_file_ch = Channel.from(params.R1)
  rev_file_ch = Channel.from(params.R2)
  
  qc_ch = QC(fw_file_ch, rev_file_ch)

  assembly_ch = assembly(qc_ch.qc_R1, qc_ch.qc_R2, params.out_fname)

  gene_file_ch = gene_prediction(assembly_ch.contigs)

  arg_ch = align_full_length_genes_with_ARG_DB(gene_file_ch.predicted_prots)

  mobileOG_ch = align_full_length_genes_with_mobileOG_DB(gene_file_ch.predicted_prots)

  metacompare(assembly_ch.contigs, gene_file_ch.predicted_genes, params.out_fname)

  output_ARGs_and_mobility(gene_file_ch.predicted_prots, arg_ch.aln_ARG_tsv, mobileOG_ch.aln_mobileOG_tsv)
  
}

