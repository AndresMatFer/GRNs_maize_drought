#!/bin/bash

#$ -l h_vmem=10G
#$ -pe serial 1


# load modules
module load bedtools


# search for the direcotries with the motif match files
cb_directory="/ngsprojects/grassgrns/results/motif_mapping/CB_mappings_filtered"
fimo_directory="/ngsprojects/grassgrns/results/motif_mapping/FIMO_mappings_filtered"


# create the new directories to store the created files
mkdir -p /ngsprojects/grassgrns/results/motif_mapping/CB_mappings_genes/
mkdir -p /ngsprojects/grassgrns/results/motif_mapping/FIMO_mappings_genes/


# loop for motif files in CB
for bed_file in $cb_directory/*
do
  filename=$(basename -- "$bed_file")
  filename="${filename%%.*}_genes.bed.gz"

  zcat $bed_file | \
  sortBed | \
  bedtools closest -a stdin -b /ngsprojects/grassgrns/data_archive/z_mays/genome/Zea_mays_genes_sorted.bed | \
  awk '{$4 = $9} 1' | \
  cut -f 1,2,3,4,5 -d ' ' | \
  sort -n -r -k5 | \
  gzip  > /ngsprojects/grassgrns/results/motif_mapping/CB_mappings_genes/$filename
  
  echo "New bed file in: /ngsprojects/grassgrns/results/motif_mapping/CB_mappings_genes/$filename"
done


# loop for motif files in FIMO
for bed_file in $fimo_directory/*
do
  filename=$(basename -- "$bed_file")
  filename="${filename%.*}_genes.bed"

  cat $bed_file | \
  sortBed | \
  bedtools closest -a stdin -b /ngsprojects/grassgrns/data_archive/z_mays/genome/Zea_mays_genes_sorted.bed | \
  awk '{$4 = $12} 1' | \
  cut -f 1,2,3,4,5 -d ' ' | \
  sort -n -r -k5  > /ngsprojects/grassgrns/results/motif_mapping/FIMO_mappings_genes/$filename

  echo "New bed file in: /ngsprojects/grassgrns/results/motif_mapping/FIMO_mappings_genes/$filename"
done

