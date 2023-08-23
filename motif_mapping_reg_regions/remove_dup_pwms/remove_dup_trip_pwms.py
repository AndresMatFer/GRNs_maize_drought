#!/usr/bin/env python3

### DESCRIPTION: REMOVE DUPLICATED AND TRIPLICATED MOTIFS
### YOU HAVE TO PROVIDE THE ARGUMENTS:
### USED FILES:
    #
### ALL FILES CAN BE FOUND IN:

import subprocess
import sys

# give the path of the Ncore.tab file
ncor_file = sys.argv[1]
# gets the path of the rest of files (should be stored in the same path as the ncor file)
path = '/'.join(ncor_file.split('/')[0:-1])

def get_related_elements(filename):

    # create two dictionaries to capture the relationships of both columns
    related_elements = {}

    with open(filename, "r") as fn:
        for line in fn:
            if line.startswith('#') or not line:
                continue

            motif1, motif2, Ncor, strand = line.split('\t')

            if Ncor == '1.000' and motif1 != motif2:
                related_elements.setdefault(motif1, set()).add(motif2)
                related_elements.setdefault(motif2, set()).add(motif1)

    return related_elements


def get_related_dict(related_elements):

    # create a set with aleady used keys
    used_keys = set()
    # create an output dict
    related_dict = {}

    # iterate over the keys (motifs) and its values (set with related motifs)
    for key, values in related_elements.items():

        # create a copy of values so we can modify it
        new_values = values.copy()

        # if key not in used keys we can iterate
        if key not in used_keys:
            
            used_keys.add(key)

            # iterate over the values (motifs)
            for value in values:

                # if value is a key in related elements merge those values with new values except the key
                if value in related_elements and value not in used_keys:

                    # first remove the key from the values
                    removed_key_set = related_elements[value].remove(key)

                    # merge the set of motifs after removing the key with new values, if there are any
                    if removed_key_set:
                        new_values = new_values.union(removed_key_set)
                    
                    # add these new key to values just to see
                    try:                        
                        for key in removed_key_set:
                            if key not in used_keys:
                                values.add(key)
                    except TypeError:
                        pass

                    # add the value to used keys so we don't use it in future iterations
                    used_keys.add(value)

            related_dict[key] = new_values

    return related_dict


def change_motif_ids(motif_info_file, related_dict):

    # create a list text to create the new text
    text = []

    with open(motif_info_file, "r") as mf:
        for line in mf:
            fields = line.rstrip().split()
            try:
                motif_id = fields[0]
            except:
                print(fields)

            if motif_id in related_dict:
                if isinstance(related_dict[motif_id], set):
                    related_dict[motif_id] = related_dict[motif_id].pop()
                fields[0] = related_dict[motif_id]

            text.append("\t".join(fields))

    return text


def remove_motif_TF_dups(text):
    no_dups_text = []
    exisiting_pairs = set()
    for line in text:
        pair = tuple(line.split()[:2])

        if pair not in exisiting_pairs:
            no_dups_text.append(line)
            exisiting_pairs.add(pair)

    return "\n".join(no_dups_text)


def main():
    motif_info_file = f"{path}/motifs_FINAL_TOTAL.txt"
    related_elements = get_related_elements(ncor_file)
    related_dict = get_related_dict(related_elements)
    text = change_motif_ids(motif_info_file, related_dict)
    no_dups_text = remove_motif_TF_dups(text)


    with open(f"{path}/motifs_FINAL_NO_DUP.txt", "w") as nodup:
        nodup.write(no_dups_text)

if __name__ == "__main__":
    main()

# import pandas as pd
# rsat_out = pd.read_csv("../matrices_comparison/rsat_output_dups_check/rsat_output_dups_check.tab", sep = "\t")
# rsat_out_ncor1 = rsat_out[rsat_out.Ncor == 1]
# duplicates = rsat_out_ncor1[rsat_out_ncor1.iloc[:, 0] != rsat_out_ncor1.iloc[:, 1]]
# duplicates.columns = ['mot1', 'mot2', 'Ncor', 'strand']
# duplicates = duplicates.drop_duplicates(subset = ['mot1', 'mot2'])
# duplicates

# rep_dict = duplicates[['mot1', 'mot2']].set_index('mot1')['mot2'].to_dict()
# mot_tf_file = "../maize_v5_motif_tf_files/maize_v5_motif_tf_duplicates.txt"
# mot_tf_df = pd.read_csv(mot_tf_file)
# mot_tf_df_nodups = mot_tf_df.replace(rep_dict).drop_duplicates()
# mot_tf_df_nodups