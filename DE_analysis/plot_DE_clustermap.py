### CREATE HEATMAP OF EXPERIMENTS TO INVESTIGATE THE RESULTS OF THE DEcalling    
### 
### USED FILES: 
#     
### ALL FILES CAN BE FOUND IN:

import os
import math
import pathlib
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


RESULTS_PATH = pathlib.Path("/ngsprojects/grassgrns/results/")
PLOTS_PATH = pathlib.Path("/ngsprojects/grassgrns/results/plots/DEG")
DEG_PATH = RESULTS_PATH / "DEG"
DOWNLOAD_PATH = pathlib.Path("/group/transreg/anmat/Descargar/DEG_plots")


def create_up_down_files(samples, outdir):

    for sample in samples:
        sample_dir = DEG_PATH / sample
        cvs_dir = sample_dir / "csv"
        contrasts_file = sample_dir / "good_contrasts_new.txt"
        with open(contrasts_file, 'r') as cf:
            contrasts = cf.read().splitlines()

        deg_files = [contrast.split("_summary.csv")[0] + '_filtered.csv' for contrast in contrasts]

        for file in deg_files:

            csv_file = cvs_dir / file
            df = pd.read_csv(csv_file, index_col=0)

            down_degs = set(df.loc[df["status_name"]=="DOWN"].index)
            up_degs = set(df.loc[df["status_name"]=="UP"].index)
            total_degs = down_degs.union(up_degs)

            # get the filenames
            upfilename = "_".join([".".join([sample] + file.split('_')[1:-3])] + ['down_degs_filtered.csv'])
            downfilename = "_".join([".".join([sample] + file.split('_')[1:-3])] + ['up_degs_filtered.csv'])
            totalfilename = "_".join([".".join([sample] + file.split('_')[1:-3])] + ['total_degs_filtered.csv'])

             # write the files
            down_file = outdir / upfilename
            up_file = outdir / downfilename
            degs_file = outdir / totalfilename

            with open(down_file, "w") as df:
                df.write('\n'.join(deg for deg in down_degs))
            
            with open(up_file, "w") as uf:
                uf.write('\n'.join(deg for deg in up_degs))

            with open(degs_file, "w") as f:
                f.write('\n'.join(deg for deg in total_degs))


# FUNCTION TO GET FEATURES
def get_features(file):
    features = [feature for feature in [file.split('_')[0].split('-')[0].split('.')[i] for i in range(7)]]
    treatment = '_'.join([feature for feature in [features[i] for i in [3,4]]])
    features = features[0:3] + [treatment.lstrip('_')] + features[5:7]
    
    features = '_'.join([feature for feature in features if feature])

    return features


# FUNCTION TO CALCULATE THE JACCARD INDEX OF TWO SETS
def calculate_jaccard_index(set1, set2):
    # get number of genes present in both sets (shared)
    nominator = len(set1.intersection(set2))
    # get total number of genes from the union of both genes
    denominator = len(set1.union(set2))
    # calculate jaccard index
    try:
        jaccard_index = nominator/denominator
    # in case the denominator is 0, the jaccard index is also 0
    except ZeroDivisionError:
        jaccard_index = 0
    
    return jaccard_index


def get_jaccard_and_features(filesdir, type):
    print(type)

    # get the list to be used as columns and rows
    col_rows = []

    # get the dictionary with features and the jaccard index
    jaccard_dict = {}

    # used combinations
    used_combs = set()

    # get all files
    files = os.listdir(filesdir)

    # iterate over every file, get features, degs and jaccard index with the rest of files
    for file1 in files:
        
        if type in file1:
            
            # get degs
            with open(filesdir / file1, 'r') as f:
                degs1 = set(f.read().splitlines())
            
            # if there are no degs skipt the dataset
            if not degs1:
                continue

            # get features of the file
            features1 = get_features(file1)
            # add them to the list of columns and rows
            col_rows.append(features1)
            
            # now get jaccard index with the rest of files
            for file2 in files:
                if type in file2:
                    
                    features2 = get_features(file2)

                    # if the combination is already used skip
                    if frozenset([features1, features2]) in used_combs:
                        continue
                    
                    # get degs
                    with open(filesdir / file2, 'r') as f:
                        degs2 = set(f.read().splitlines())
                    
                    # if there are no degs skipt the dataset
                    if not degs2:
                        continue

                    # calculate jaccard index
                    jaccard_index = calculate_jaccard_index(degs1, degs2)

                    # add the jaccard index to the jaccard dict
                    jaccard_dict[frozenset([features1, features2])] = jaccard_index

                    # add the combination to used combinatons so we don't calculate the jaccard index again
                    used_combs.add(frozenset([features1, features2]))
    
    return col_rows, jaccard_dict


# PLOT CLUSTERMAP
def plot_clustermap(col_rows, jaccard_dict, deg_type, out_path=PLOTS_PATH, download_path=DOWNLOAD_PATH):
    # variable for later use in clustermap
    vmax = 0
    # create data frames for clustermap
    df = pd.DataFrame()
    
    columns = col_rows.copy()
    rows = col_rows.copy()

    # get tissue for the color in the clustermap
    tissue = [column.split("_")[1] for column in columns]
    
    tissue_unique = []
    for column in columns:
        if column.split("_")[1] not in tissue_unique:
            tissue_unique.append(column.split("_")[1])
    tissue_unique.sort()

    for column in columns:
        for row in rows:
            jaccard_val =  jaccard_dict[frozenset([column, row])]
            df.loc[row, column] = jaccard_val
            if jaccard_val > vmax and jaccard_val != 1:
                vmax = math.floor(jaccard_val*10)/10 # to make the scale always a bit bigger

    vmax += 0.1

    # plot the clustermaps
    plt.figure(figsize=(60, 50))
    lut = dict(zip(tissue_unique, ["orange", "b", "y", "g",  "c", "purple", "r", "m", "k"]))    
    row_colors = [lut[s] if any(s in label for label in df.index) else None for s in tissue]
    sns.set(font_scale=2.4)
    clustermap = sns.clustermap(df, vmax=0.5, row_colors=row_colors, figsize=(60,50), robust=True)
    # clustermap.set_xticklabels(clustermap.get_xticklabels(), rotation=90)  # Rotate the x labels
    # cbar = clustermap.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=18)
    plt.xticks(fontsize=30)  
    plt.yticks(fontsize=30)
    # plt.title(f"Clustermap {deg_type} DEGs with Jaccard index", fontdict={'fontsize': 80, 'weight': 'bold'})
    plt.tight_layout()
    plt.savefig(out_path / f"cluster_map_{deg_type}_quantiles.pdf", format="pdf")
    plt.savefig(download_path / f"cluster_map_{deg_type}_quantiles.pdf", format="pdf")


def main():
    # create a new dir in results/DEG_up_down_csv
    up_down_csv = RESULTS_PATH / "DEG_up_down_csv" / "all_csv"
    up_down_csv.mkdir(exist_ok=True)

    # for every file in good_contrasts create an up, down and total file in the new directory
    samples = os.listdir(DEG_PATH)

    # write files
    # create_up_down_files(samples, up_down_csv)

    # iterate three times, if they have up, down or total in the name of the file
    # for every file save the info sample_tissue_devstage_treatment_genotype_time in features
    # then for every file read every other file and get jaccard index, store it in a dict such as: {frozenset(feat1, feat2):jaccard}
    # every feature should be unique, so we can store both features in a set and if they are already present we don't get the jaccard_idx
    # also for columns and rows every feature we store it in a list to be used later to create the cluster_map 
    up_col_rows, up_jaccard_dict = get_jaccard_and_features(filesdir=up_down_csv, type = "up_degs")
    down_col_rows, down_jaccard_dict = get_jaccard_and_features(filesdir=up_down_csv, type = "down_degs")
    total_col_rows, total_jaccard_dict = get_jaccard_and_features(filesdir=up_down_csv, type = "total_degs")

    # plot clustermaps
    # plot_clustermap(up_col_rows, up_jaccard_dict, deg_type="UP", out_path=PLOTS_PATH)
    # plot_clustermap(down_col_rows, down_jaccard_dict, deg_type="DOWN", out_path=PLOTS_PATH)
    plot_clustermap(total_col_rows, total_jaccard_dict, deg_type="TOTAL", out_path=PLOTS_PATH)


if __name__ == '__main__':
    main()