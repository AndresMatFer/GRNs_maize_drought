#!/bin/bash

feat_files_path="/ngsprojects/grassgrns/results/ensemble_motif_mapppins_nico_coords_TG/feature_files"
degs_file=$1 	# file with two columns (set_id	  deg)					("/ngsprojects/grassgrns/results/DEG_selected_sets/degs_in_leaf_cluster.txt")
out_dir=$2	# directory where you want the enrichment output files to be created 	("/ngsprojects/grassgrns/results/motif_enrichment/leaf_enrichment")
key_word=$3	# something like root or leaf 						("leaf")

# enrichment command
for file in $feat_files_path/* 
do
	# get the names of the new files
	filename=$(basename "$file")
	new_file="${filename%_feature_file_no_dup.txt}"
	
	# perform enrichment
	/group/transreg/Tools/dreec_enrichment/enricherv2.4 -p \
	 "$file" \
	 "$degs_file" \
	 > "${out_dir}/${new_file}_${key_word}.txt"
	
done	


#/group/transreg/Tools/dreec_enrichment/enricherv2.4 \
#$file \ /validation_for_enricher.txt > test_enricher_output
