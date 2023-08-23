### GENERATE THREE NETWORK FILES TO BE USED IN GO ENRICHMENT (one for UP, one for DOWN and one for TOTAL)
### 
### USED FILES: maize_15up_2.5down_acrs_metadata.txt
#     
### ALL FILES CAN BE FOUND IN:

import os
import sys
import pathlib
import pandas as pd

OUT_DIR =   pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment/network_files_for_go/degs_15up_2.5down_acrs")


def write_network_file(df:pd.DataFrame, filename:None):
    tf_tgs_dict = {}
    for _, row in df.iterrows():

        tf = row["gene_id"]
        tgs = set(row["hits"].split(","))

        tf_tgs_dict.setdefault(tf, set()).update(set(tgs))
        

    lines = []
    for tf, tgs in tf_tgs_dict.items():
        for tg in tgs:
            lines.append(f"{tf}\t{tg}")

    with open(OUT_DIR / f"network_file_{filename}_15up_2.5down_acrs.txt", "w") as out:
        out.write('\n'.join(lines))


def get_full_df_network_file(df:pd.DataFrame, curated:None):

    df = df[df["curated"] == curated]

    up_set = set(df[df["set_id"]=="UP"]["gene_id"])
    down_set = set(df[df["set_id"]=="DOWN"]["gene_id"])

    inter_set = up_set.intersection(down_set)
    # just keep the TFs unique for up and down
    up_set = up_set.difference(inter_set)
    down_set = down_set.difference(inter_set)

    print(f"up {len(up_set)}")
    print(f"down {len(down_set)}")
    print(f"inter {len(inter_set)}")

    up_df = df[df['gene_id'].isin(up_set)]
    down_df = df[df['gene_id'].isin(down_set)]
    inter_df = df[df['gene_id'].isin(inter_set)]

    name_curated = {0:"FPs", 1:"TPs"}

    write_network_file(up_df, filename=f"{name_curated[curated]}_UP")
    write_network_file(down_df, filename=f"{name_curated[curated]}_DOWN")
    write_network_file(inter_df, filename=f"{name_curated[curated]}_INTER")



def main():
    # file = sys.argv[1]
    file = "/ngsprojects/grassgrns/results/motif_enrichment/enrichment_metadata/degs_in_10_contr/maize_15up_2.5down_acrs_degs_in_10_samples_metadata.txt"
    df = pd.read_csv(file, sep="\t")
    # first we sort the df by rank
    df = df.sort_values("rank", ascending=True)
    # select only the rows with rank <= 100
    df = df[df['rank'] <= 100]

    get_full_df_network_file(df, curated=0)
    get_full_df_network_file(df, curated=1)


    
if __name__ == "__main__":
    main()