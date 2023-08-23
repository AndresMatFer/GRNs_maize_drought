### PLOT BARCHART FOR TOP GO TERMS 
### USED FILES: files from the output of the go terms enricher
#     
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/results/go_enrichment

import os
import sys
import pathlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

GO_PATH = pathlib.Path("/ngsprojects/grassgrns/results/go_enrichment/degs_15up_2.5down_acrs")
METADATA_PATH = pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment/enrichment_metadata")
OUT_PATH = pathlib.Path("/ngsprojects/grassgrns/results/plots/go_plots/degs_15up_2.5down_acrs/curated")
DOWNLOAD_PATH = pathlib.Path("/group/transreg/anmat/Descargar/go_enrichment_output/plots/degs_15up_2.5down_acrs/curated")


def get_athaliana_orthologs_and_description(metadata_file):
    meta_df = pd.read_csv(metadata_file, sep='\t', header=0)
    zm_at_dict = {}
    for _,row in meta_df.iterrows():
        try:
            name_fig = "(" + " - ".join(row["gene_name_figs"].rstrip(")").split(" (")) + ")"
        except AttributeError:
            name_fig = " "
            
        description = row["description"]
        curated = row["curated"]
        if curated:
            zm_at_dict[row["gene_id"]] = f"{description} CURATED {name_fig}"
        else:
            zm_at_dict[row["gene_id"]] = f"{description} {name_fig}"

    # zm_at_dict = {row["gene_id"]:row["ath_alias"] for _,row in meta_df.iterrows()}

    return zm_at_dict


def get_meta_15up_2down_acrs_files(metadata_path=METADATA_PATH):
    meta_files = [
        metadata_path / path / file
        for path in os.listdir(metadata_path)
        for file in os.listdir(metadata_path / path)
        if "5up_1down_acrs" in file
    ]

    return meta_files

    
def get_tfs(file):

    df = pd.read_csv(file, sep='\t')

    # first we sort the df by rank
    df = df.sort_values("rank", ascending=True)

    # select only the rows with rank <= 100
    df_100 = df[df['rank'] <= 100]
    # and select only the rows with TPs (curated = 1)
    df_tps = df[df['curated'] == 1]

    dict_100 = {}
    # then we get the TFs for the rank100 motifs
    dict_100["UP"] = set(df_100[df_100["set_id"] == "UP"]["gene_id"])
    dict_100["DOWN"] = set(df_100[df_100["set_id"] == "DOWN"]["gene_id"])
    dict_100["TOTAL"] = set(df_100[df_100["set_id"] == "TOTAL"]["gene_id"])
    
    dict_tps = {}
    # and for the curated motifs
    dict_tps["UP"] = set(df_tps[df_tps["set_id"] == "UP"]["gene_id"])
    dict_tps["DOWN"] = set(df_tps[df_tps["set_id"] == "DOWN"]["gene_id"])
    dict_tps["TOTAL"] = set(df_tps[df_tps["set_id"] == "TOTAL"]["gene_id"])

    return dict_100, dict_tps


def create_heatmap(filename, tfs_dict:dict, tfs_in:str, zm_at_dict:dict):
    degtype = filename.split("_")[2]
    new_name = filename.split("network_file_")[1].split(".txt")[0]
    title = filename.split("network_file_")[1].split("_15up_2.5down_acrs_GO_enrichment.txt")[0]
    title = " ".join(title.split("_"))

    filepath = GO_PATH / filename

    df = pd.read_csv(filepath, sep='\t', header=None)
    
    # we select the tfs of the df that are in the dictionary
    selected_tfs = set(df[0]).intersection(tfs_dict[degtype])

    # we want a 2d object that is the selected tfs rows, the go terms columns
    selected_df = df[df[0].isin(selected_tfs)]
    selected_go_terms = set(selected_df[9])

    heat_map_df = pd.DataFrame(0, index=selected_go_terms, columns=selected_tfs)

    # put the q-values and if there is no association put 0
    for row in selected_df.itertuples(index=False):
        tf = row._0
        go = row._9
        p_val = row._2
        q_val = row._3
        enr_fold = row._4

        # pi_value (Xiao et al. 2014) --> pi_value = enr_fold * -log_10(p-val)
        pi_val = enr_fold*-np.log10(p_val)

        heat_map_df.loc[go, tf] = -np.log10(q_val)    # option for -log10(qval)
        # heat_map_df.loc[go, tf] = pi_val              # option for pi-value
    
    assert tfs_in in ["rank100", "TPs"], "tfs_in must be a 'rank100' or 'TPs'"

    heat_map_df.columns = [zm_at_dict[tf] for tf in selected_tfs]

    # I decided to take out go terms that are not associated with at least 10% of TFs
    percentages = []
    for _,row in heat_map_df.iterrows():
        num_zeros = (row == 0).sum() 
        percentage_zeros = num_zeros / len(row) 
        percentages.append(percentage_zeros)
    
    # Filter the DataFrame to exclude rows with percentage of zeros >= 20%
    # print(percentages)
    percentages = np.array(percentages)
    heat_map_df = heat_map_df[percentages < 0.80]

    # make the plot
    sns.set(style="whitegrid", font_scale=2.0)
    plt.figure(figsize=(60,50))
    # sns.set(font_scale=0.8)
    clustermap = sns.clustermap(heat_map_df, figsize=(60,50), yticklabels=True ,robust=True)
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right")
    # plt.xticks(rotation = 45,ha='right')
    plt.title(title, fontdict={"fontsize":18}, loc="center")
    # plt.yticks(fontsize=40)
    # plt.yticks(fontsize=40)
    # clustermap.ax_heatmap.xaxis.set_tick_params(labelsize=4)
    # clustermap.ax_heatmap.yaxis.set_tick_params(labelsize=4)
    plt.tight_layout()
    plt.savefig(OUT_PATH / f"clustermap_{new_name}_tfs_in_{tfs_in}_curated.pdf", format="pdf")
    plt.savefig(DOWNLOAD_PATH / f"clustermap_{new_name}_tfs_in_{tfs_in}_curated.pdf", format="pdf")
    plt.close()



def main():
    print("working")
    meta_files = get_meta_15up_2down_acrs_files(metadata_path=METADATA_PATH)  
    go_files = os.listdir(GO_PATH)  
    print(go_files)
    for file in meta_files:
        study = str(file).split('/')[6].split("_metadata")[0]
        print(study)
        zm_at_dict = get_athaliana_orthologs_and_description(file)
        dict_100, dict_tps = get_tfs(file)

        for go_file in go_files:
            if study.split("_")[0] not in go_file:
                continue
            # if "_degs_" not in study:
            #     continue
            print(file)
            print(go_file)
            create_heatmap(go_file, dict_100, tfs_in="rank100", zm_at_dict=zm_at_dict)
            create_heatmap(go_file, dict_tps, tfs_in="TPs", zm_at_dict=zm_at_dict)


    # # get tfs dictionaries
    # meta_files = get_meta_1up1down_acrs_files(metadata_path=METADATA_PATH)    
    # file = meta_files[0]
    # print(file)
    # dict_100, dict_tps = get_tfs(file)
    
    # # 

    # go_files = os.listdir(GO_PATH)
    # filename = go_files[-1]
    # print(filename)
    # create_heatmap(filename, dict_100)

if __name__ == "__main__":
    main()