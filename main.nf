#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  if(!params.fastaSubsetSize) {
    error "Missing params.fastaSubsetSize"
  }

  if(params.inputFilePath) {
    seqs = channel.fromPath( params.inputFilePath )
             .splitFasta( by:params.fastaSubsetSize, file:true  )
  }
  else {
    error "Missing params.inputFilePath"
  }

  tmhmm(seqs)
  gff = tmhmm2gff(tmhmm.out)
  indexResults(gff.collectFile(), params.outputFileName)
}

process tmhmm {

input:
    path subsetFasta

  output:
    path 'tmhmm_results'

  script:
  """
  tmhmm -short $subsetFasta >tmhmm_results
  """
}


process tmhmm2gff {
  container 'bioperl/bioperl:stable'

  input:
    path subset

  output:
    path 'tmhmm_subset.gff'

  script:
  """
  tmhmm2gff.pl --inputFile $subset \
                 --outputFile tmhmm_subset.gff
  """

}


process indexResults {
  container 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gz'
    path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n -k5,5nr $gff > ${outputFileName}
  bgzip ${outputFileName}
  tabix -p gff ${outputFileName}.gz
  """
}
