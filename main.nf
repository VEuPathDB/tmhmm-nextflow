#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.inputFilePath) {
  seqs = Channel.fromPath( params.inputFilePath )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.inputFilePath")
}

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  tmhmm(seqs)
    | collectFile(name: params.outputFileName, storeDir: params.outputDir)
}

process tmhmm {

input:
    path subsetFasta

  output:
    path 'tmhmm_results'

  script:
  """
  tmhmm $task.ext.args $subsetFasta >tmhmm_results
  """
}
