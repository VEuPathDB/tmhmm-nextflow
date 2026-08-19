# tmhmm-nextflow

Nextflow pipeline that predicts transmembrane helices in protein sequences using TMHMM and publishes the results as an indexed GFF3 file.

## Overview

This pipeline runs [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/) over a set of input protein sequences to identify transmembrane helix regions and their topology (inside/outside orientation relative to the membrane). It is used within VEuPathDB's genome annotation workflows to generate transmembrane domain evidence for predicted proteins, which feeds into downstream protein feature annotation.

The workflow splits the input FASTA file into subsets, runs TMHMM on each subset in parallel, converts the raw TMHMM output into GFF3 features (including `TMhelix`, `inside`, and `outside` subfeatures per protein), then merges, sorts, `bgzip`-compresses, and `tabix`-indexes the combined result.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker or Singularity, depending on the execution profile in the runtime configuration

## Usage

```
nextflow run VEuPathDB/tmhmm-nextflow \
  -r main \
  -resume \
  --inputFilePath /path/to/proteins.fa \
  --fastaSubsetSize 500 \
  --outputFileName tmhmm.out \
  --outputDir /path/to/output \
  -C <config>
```

The pipeline has a single, unnamed entry point (no `-entry` flag is needed).

## Key parameters

| Parameter | Description |
| --- | --- |
| `inputFilePath` | Path to the input protein FASTA file. Required. |
| `fastaSubsetSize` | Number of sequences per FASTA subset chunk processed in each parallel TMHMM job. Required. |
| `outputFileName` | Name of the merged, sorted GFF3 output file (default `tmhmm.out`). |
| `outputDir` | Directory the final compressed/indexed output is published to (default `output` under the launch directory). |

Container images and process-level resources (executor, queue, memory) are supplied via the runtime configuration passed with `-C`, allowing the same pipeline definition to run under Docker, Singularity, or an LSF cluster profile (see `conf/docker.config`, `conf/singularity.config`, and `conf/lsf.config`).

## Output

For each run, the pipeline publishes to `outputDir`:

- `<outputFileName>.gz` — a sorted GFF3 file with one feature per input protein describing predicted transmembrane helices and their inside/outside membrane topology (attributes include `ExpectedAA`, `First60`, and `PredictedHelices` counts from TMHMM's short-format output).
- `<outputFileName>.gz.tbi` — a Tabix index for the compressed GFF3 file.
