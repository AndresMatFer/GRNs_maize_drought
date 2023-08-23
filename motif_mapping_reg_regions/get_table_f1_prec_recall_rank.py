### CALCULATE F1, PRECISSION, RECALL OF THE ENRICHMENT FILES AND THE CURATED VALIDATION FILE
### 
### USED FILES: 
#     
### ALL FILES CAN BE FOUND IN:

import os
import pathlib
import pandas as pd

ENRICH_PATH = pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment")
METADATA_PATH = pathlib.Path("/ngsprojects/grassgrns/results/motif_enrichment/enrichment_metadata")
DOWNLOAD_PATH = pathlib.Path("/group/transreg/anmat/Descargar/enrichment_output")


def calc_precision_recall_f1(curated_set, motifs_set):
    tp = motifs_set.intersection(curated_set) # motifs in both sets
    fp = motifs_set.difference(curated_set)   # motifs present in degs_in_10_samples that are not in curated
    fn = curated_set.difference(motifs_set)   # motifs present in curated that are not in degs_in_10_samples

    try:
        precision = round(len(tp) / len(tp.union(fp)), 4)
    except ZeroDivisionError:
        precision = 0
    try:
        recall    = round(len(tp) / len(tp.union(fn)), 4)
    except ZeroDivisionError:
        recall    = 0
    try:
        f1_score  = round((2 * precision * recall) / (precision + recall), 4)
    except ZeroDivisionError:
        f1_score   = 0
    num_tp = len(tp)

    return precision, recall, f1_score, num_tp


def get_rank(filename:str, meta_subpath):
    
    name = filename#.split(".txt")[0] + "_metadata.txt"
    meta_file = METADATA_PATH / meta_subpath / name 

    df = pd.read_csv(meta_file, sep='\t')

    ranks = dict()

    for _,row in df.iterrows():
        if row["curated"] == 1 and row["set_id"] == "TOTAL":
            ranks.setdefault(row["gene_id"], set()).add(row["rank"]) # we are interested in the minimum TF rank associated
    
    ranks_list = [min(ranks[gene_id]) for gene_id in ranks.keys()]

    return round(sum(ranks_list) / len(ranks_list), 4)  


def main():
    TABLES_OUTDIR = ENRICH_PATH / "promoter_metric_tables"
    curated_file = "/ngsprojects/grassgrns/results/ensemble_motif_mapppins_nico_coords_TG/curated_feature_file_for_validation.txt"

    df_lines = []
    
    df = pd.read_csv(curated_file, sep='\t', header=None)
    curated_set = set(df[0])

    subpath = "root_degs_5_contrasts_metadata" # "degs_in_10_samples" "leaf_metadata" "leaf_degs_8_contrasts_metadata" "root_metadata" "root_degs_5_contrasts_metadata"

    for filename in os.listdir(METADATA_PATH / subpath):
        file = METADATA_PATH / subpath / filename
        df = pd.read_csv(file, sep='\t')
        # df = pd.read_csv(file, sep='\t', comment='#', header=None)
        
        # df.columns = ["set_id","ftr_id","p-val","q-val","enr_fold","set_size","ftr_size","n_hits","hits"]
        # only total
        df = df[df["set_id"] == "TOTAL"]
        # CHANGE THE Q-VALUE AND THE ENRICHEMNT FOLD THRESHOLD TO BE MORE RESTRICTIVE
        q_val_threshold = 0.01    # 0.01, 0.001, 0.0001
        enr_fold_threshold = 1.2    # based on the enr_fold distribution plots
        df = df[df["q-val"] < q_val_threshold]
        df = df[df["enr_fold"] > enr_fold_threshold]


        motifs_set = set(df["ftr_id"])
        
        precision, recall, f1_score, tp_set = calc_precision_recall_f1(curated_set, motifs_set) 
        rank = get_rank(filename, meta_subpath=subpath)
        tp_total = len(curated_set)

        filename = "_".join(filename.split("_")[1:-1])
        num_motifs = len(motifs_set) # change to total_set
        df_lines.append([filename,precision,recall,f1_score,rank,num_motifs,tp_set,tp_total])
    
    df = pd.DataFrame(df_lines)
    df.columns = ["promoter","precision","recall","f1 score","average min rank","number of motifs","TPs in set","TPs total"]
    df = df.sort_values("f1 score", ascending=False)

    key_word = subpath.split("_metadata")[0] # can be leaf, root or degs_in_10_samples

    df.to_csv(TABLES_OUTDIR / f"{key_word}_promoter_metrics_q_val_{q_val_threshold}_enr_fold_{enr_fold_threshold}.tsv", sep='\t', index=False)

    with pd.ExcelWriter(DOWNLOAD_PATH / f"{key_word}_promoter_metrics_q_val_{q_val_threshold}_enr_fold_{enr_fold_threshold}.xls") as writer:
        df.to_excel(writer, index = False)


if __name__ == "__main__":
    main()