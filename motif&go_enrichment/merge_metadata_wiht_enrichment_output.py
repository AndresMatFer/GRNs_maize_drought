### MERGE METADATA INFO FILE OF NICO WITH THE OUTPUT OF THE ENRICHER
###
### USED FILES: we need to use the motifs_FINAL_NO_DUP.txt
#
### ALL FILES CAN BE FOUND IN:

import os
import pathlib
import numpy
import pandas as pd

ENRICH_PATH = pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment/no_metadata")
OUT_PATH = pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment/enrichment_metadata/root_metadata_degs_5_contrasts")
DOWNLOAD_PATH = pathlib.Path("/group/transreg/anmat/Descargar/enrichment_output")


def assign_metadata(enrich_file, metadata_file, motif_TF_file, curated_file):

    # create dataframe from the motif enrichment file
    df = pd.read_csv(enrich_file, sep='\t', comment='#', header=None)
    column_names = ["set_id", "ftr_id", "p-val", "q-val", "enr_fold", "set_size", "ftr_size", "n_hits","hits"]
    df.columns = column_names

    # ranking by pi_value (Xiao et al. 2014) --> pi_value = enr_fold * -log_10(p-val)
    df["pi-val"] = df["enr_fold"]*(-numpy.log10(df["p-val"]))
    new_columns_order = ["set_id", "ftr_id", "p-val", "q-val", "enr_fold", "pi-val", "set_size", "ftr_size", "n_hits", "hits"]
    df = df.reindex(columns=new_columns_order)
    df["rank"] = df.groupby("set_id")["pi-val"].rank(method="dense", ascending=False)

    ## TO ADD THE METADATA (we need the motif-TF file and the metadata file)
    # get metadata dataframe
    meta = pd.read_csv(metadata_file, sep=",")    

    # get dictionary with motif_id:set_of_TFs
    motif_tf_dict = {}
    with open(motif_TF_file, 'r') as mtf:
        for line in mtf:
            if line.startswith("motif_id"):
                continue  # skip 1st line 
            motif = line.rstrip().split()[0].split('.')[0] # some motifs en with .1
            tf = line.rstrip().split()[1]
            motif_tf_dict.setdefault(motif, set()).add(tf)

    # get list of motif-TF from the curated file (experimental)
    with open(curated_file, 'r') as c:
        curated_motif_tf = [line.rstrip().split() for line in c.readlines()]

    
    # create new df with metadata and check if curated
    new_df = []
    for _, row in df.iterrows():
        motif = row["ftr_id"]

        if motif in motif_tf_dict:
            new_rows = []
            gene_set = motif_tf_dict[motif]

            for gene in gene_set:
                for _,metarow in meta.loc[meta["gene_id"]==gene].iterrows():
                    new_row = row.copy()
                    new_row = pd.concat([new_row, metarow], axis=0)

                    if [motif, gene] in curated_motif_tf:
                        new_row["curated"] = 1
                    else:
                        new_row["curated"] = 0

                    new_rows.append(new_row)
                
            new_df += new_rows

        else:
            new_columns = list(df.columns) + list(meta.columns)
            row = row.reindex(index=new_columns)
            row["curated"] = 0
            new_df.append(row)


    # Create a new DataFrame from the list of rows
    new_df = pd.DataFrame(new_df)

    # Reset the index of the new DataFrame
    new_df.reset_index(drop=True, inplace=True)

    # add hits at the end so they don't disturb
    final_columns_order = list(df.columns)[:-1] + list(meta.columns) + ["hits"]
    df = df.reindex(columns=final_columns_order)

    return new_df


def main():
    motif_TF_file = "/ngsprojects/grassgrns/data_archive/z_mays/motif_info/motifs_FINAL_NO_DUP.txt"
    metadata_file = "/group/transreg/niman/miniac_repo/mini_ac_internal_inputs/internal_inputs/zma_v5/maize_v5_gene_metadata_file.txt"
    curated_file = "/ngsprojects/grassgrns/results/ensemble_motif_mapppins_nico_coords_TG/curated_feature_file_for_validation.txt"

    enrich_folder = "root_enrichment_degs_5_contrasts"   # "root_enrichment" "leaf_enrichment" "leaf_enrichment_degs_8_contrasts" "root_enrichment_degs_5_contrasts" "degs_in_10_samples" 
    filenames = os.listdir(ENRICH_PATH / enrich_folder)
    print(filenames)
    down_folder = "tables_root_degs_5_contrasts" # "tables_root" "tables_leaf" "tables_degs_in_10_samples" "tables_leaf_degs_8_contrasts" "tables_root_degs_5_contrasts"
    download_path = DOWNLOAD_PATH / down_folder

    for filename in filenames:
        enrich_file = ENRICH_PATH / enrich_folder / filename
        outdf = assign_metadata(enrich_file, metadata_file, motif_TF_file, curated_file)

        outfile = filename.split(".txt")[0] + "_metadata.txt"
        excelfile = filename.split(".txt")[0] + "_metadata.xls"

        outdf.to_csv(OUT_PATH / outfile, sep='\t', index=False)

        with pd.ExcelWriter(download_path / excelfile) as writer:
            outdf.to_excel(writer, index = False)


if __name__ == "__main__":
    main()