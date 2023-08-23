#!/usr/bin/env python3

### GET PWMs IN CLUSTER BUSTER FORMAT FROM PWMs IN MEME FORMAT STORED IN JSON FORMAT
### YOU HAVE TO PROVIDE THE ARGUMENTS: motif_pwm_path
### USED FILES: 
    #  motifs_FINAL_TOTAL --> For both maize and wheat constructed previously with fields: motif_id, tf_gene_id, species, pwm_json
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/

import sys
import json
import pandas as pd
import pathlib


MOTIF_FILE = sys.argv[1]


# CONVERT A PWM FROM MEME FORMAT (stored in json format) TO CLUSTER BUSTER
def get_cluster_buster(meme_pwm_dict):
    meme_pwd_json = pd.read_json(json.dumps(meme_pwm_dict))

    cb_pwd = ''
    for index, row in meme_pwd_json.iterrows():
        a = str(int(row['A']*100))
        c = str(int(row['C']*100))
        g = str(int(row['G']*100))
        t = str(int(row['T']*100))
        new_row = '\t'.join([a, c, g, t])

        cb_pwd += new_row + '\n'
    
    return cb_pwd.rstrip()
    
    
# WRITE FILE/FILES IN CLUSTER BUSTER FORMAT PWMs
def write_cluster_buster(motif_file_path, outfolder):
    
    FILES_DIR = pathlib.Path('/'.join(motif_file_path.split("/")[:-1]))

    with open(motif_file_path, "r") as mf:
        text = ''
        next(mf)

        if outfolder:
            outfolder = FILES_DIR / outfolder
            outfolder.mkdir(parents=True, exist_ok=True)

        motifs_in_all = set()
        for line in mf:
            fields = line.rstrip().split()

            motif_id = fields[0]
            # gene_id = fields[1]
            # species = fields[2]
            pwm_dict = eval(fields[3])

            cb_pwm = get_cluster_buster(pwm_dict)
            
            text += f">{motif_id}\n{cb_pwm}\n"
            motifs_in_all.add(motif_id)

            with open(outfolder / f"{motif_id}.cb", "w") as cb:
                cb.write(f">{motif_id}\n{cb_pwm}\n")


def main():
    outfolder = "cb_pwms"

    write_cluster_buster(
        motif_file_path=MOTIF_FILE,
        outfolder=outfolder
    )

if __name__ == '__main__':
    main()