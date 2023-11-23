# giraffe-svs-nf

Genotyping structural variants with pangenomes

## Introduction

This is a nextflow pipeline for genotyping structural variants from short reads
using a pangenome reference and various tools in the [VG toolkit][vg].

It requires the following inputs:
* giraffe indices for a pangenome reference graph, including gbz, min, and dist
* paired-end short reads for a set of sample you want to genotype
* the reference to genotype against, in fasta format --- must be the same one
  chosen as the reference when building the pangenome
* a sample sheet in csv format giving sample IDs and absolute paths to fastqs
  for each sample to be genotyped (see [Preparing data](#preparing-data))

and produces a single gzipped vcf containing genotype calls for all of the
samples.

This pipeline is based on the process used to genotype SVs in [Sirén et al.
(2021)](siren).

## Installation

### On Lewis

If you are using the Lewis cluster and are a member of the `warrenlab-group`,
you can just install nextflow in a conda environment if you haven't already,
and then run the pipeline with the argument `-profile lewis`.

### With conda

If you have conda installed and in your bin, nextflow can create an
environment to run all tasks in if you run the pipeline with the argument
`-profile conda`.

### Requirements

If using one of the methods above, you don't need to worry about installing the
requirements yourself, but here they are in case you do:
* [VG toolkit][vg]
* [bcftools][bcf]
* [htslib][hts]

## Preparing data

You'll need to make a sample sheet in csv format with columns named
`sample_id`, `r1`, and `r2`. The sample ID should be a unique identifier for
that individual, and absolute paths to forward and reverse fastq files should
be specified in the r1 and r2 columns, respectively. If there are multiple
pairs of forward and reverse fastq file for a single sample, you can have
multiple rows for that sample.

Here is an example sample sheet:
```
sample_id,r1,r2
princess,/data/reads/princess_R1.fastq.gz,/data/reads/princess_R2.fastq.gz
junebug,/data/reads/junebug_R1.fastq.gz,/data/reads/junebug_R2.fastq.gz
oreo,/data/reads/oreo_L1_R1.fastq.gz,/data/reads/oreo_L1_R2.fastq.gz
oreo,/data/reads/oreo_L2_R1.fastq.gz,/data/reads/oreo_L2_R2.fastq.gz
freddie,/data/reads/freddie_R1.fastq.gz,/data/reads/freddie_R2.fastq.gz
```

## Running

To run the pipeline, use the following command:
```bash
nextflow run WarrenLab/giraffe-svs-nf \
     --sample-sheet samples.csv \
     --gbz-index index/Amex-pg.full.gbz \
     --dist-index index/Amex-pg.full.dist \
     --min-index index/Amex-pg.full.min \
     --ref-fasta index/surface.fa
     --ref-path surface
```
`--ref-path` is the name of the reference path in the pangenome graph, as
specified when building it.

**N.B.** You do not need to download or clone the code in this repository. The
command above tells nextflow fetch the pipeline code from this repository and
run it.

[siren]: https://www.science.org/doi/10.1126/science.abg8871

[vg]: https://github.com/vgteam/vg

[bcf]: https://github.com/samtools/bcftools

[hts]: https://github.com/samtools/htslib
