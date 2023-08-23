#!/bin/bash

#$ -l h_vmem=5G
#$ -pe serial 4

module load nextflow

dataset=$1

# requires run nf-core/rnaseq to be installed in the current directory
nextflow  run nf-core/rnaseq \
--input /scratch/grassgrns/data/z_mays_drought/datasets_batch_4/$dataset/samplesheet_$dataset.csv \
--outdir /scratch/grassgrns/analysis/$dataset \
--multiqc_title $dataset \
--pseudo_aligner salmon \
-profile singularity \
-c /scratch/grassgrns/data/nfcore_rna_seq_plaza.config \
-w /ngsprojects/grassgrns/data_archive/work_$dataset \
--salmon_quant_libtype A -r 3.9 -resume && rm -r /ngsprojects/grassgrns/data_archive/work_$dataset
