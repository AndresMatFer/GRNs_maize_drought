#!/bin/bash

#$ -l h_vmem=8G
#$ -pe serial 6

module load salmon
salmon index -i /ngsprojects/grassgrns/data_archive/z_mays/genome/zea_mays_salmon_index_decoys --threads 6  --kmerLen 31 --tmpdir /ngsprojects/grassgrns/data_archive/salmon_index/tmp --keepFixedFasta --transcripts /ngsprojects/grassgrns/data_archive/z_mays/genome/Zea_mays_GENTROME.fa --keepDuplicates -d /ngsprojects/grassgrns/data_archive/z_mays/genome/decoys.txt
