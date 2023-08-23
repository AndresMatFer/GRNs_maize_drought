#!/bin/bash

module load bedtools

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_bed_file> <keyword>"
  exit 1
fi

# Assign the input arguments to variables
input_bed_file="$1"
keyword="$2"

# The regulatory region bed files are sotored here
reg_bed_path="/ngsprojects/grassgrns/data_archive/z_mays/reg_regions/bedfiles"

# Define the output file name based on the input bed file and keyword
output_bed_file="${input_bed_file%.*}_${keyword}.bed"

# Determine the second bed file based on the keyword
case $keyword in
  acrs)
    second_bed_file="$reg_bed_path/ACRs-combined_all-organs_sorted.bed"
    ;;
  umrs)
    second_bed_file="$reg_bed_path/zea_mays_umrs_v5_sorted.bed"
    ;;
  cnss)
    second_bed_file="$reg_bed_path/zea_mays_cnss_v5_sorted.bed"
    ;;
  *)
    echo "Invalid keyword. Available options are: acrs, umrs, cnss"
    exit 1
    ;;
esac

# Perform bedtools intersect using the input and second bed files
bedtools intersect -a "$input_bed_file" -b "$second_bed_file" > "$output_bed_file"

echo "Intersection completed. Output stored in: $output_bed_file"

