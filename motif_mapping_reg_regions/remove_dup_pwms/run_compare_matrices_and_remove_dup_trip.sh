#$ -l h_vmem=20G
#$ -M andresbiotec99@gmail.com
#$ -m bes

module load rsa-tools
module load perl/x86_64/5.14.1
module load python

## input file with all PWMs in cluster-buster format (/ngsprojects/grassgrns/data_archive/z_mays/motif_info/Zea_mays_all_motifs_cb.cb)
cb_file=$1

## name of output file with PWMs in transfac format (/ngsprojects/grassgrns/data_archive/z_mays/motif_info/Zea_mays_all_motifs_transfac.txt)
transfac_file=$2

## we use transfac format because cluster-buster does not work as expected
## so first is conversion of cluster-buster PWMs from cb to transfac and then compare
convert-matrix -i $cb_file -o $transfac_file -from cluster-buster -to transfac

compare-matrices -file $transfac_file -format transfac -o /ngsprojects/grassgrns/data_archive/z_mays/motif_info/Z_mays_Ncor -lth Ncor 0

# run the remove duplicates and triplicates
python3 remove_dup_trip_pwms.py /ngsprojects/grassgrns/data_archive/z_mays/motif_info/Z_mays_Ncor.tab
