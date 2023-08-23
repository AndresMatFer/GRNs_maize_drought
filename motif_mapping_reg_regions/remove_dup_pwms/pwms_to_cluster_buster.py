#!/usr/bin/env python3

### GET PWMs IN CLUSTER BUSTER FORMAT FROM PWMs IN MEME FORMAT STORED IN JSON FORMAT
### YOU HAVE TO PROVIDE THE ARGUMENTS: motif_pwm_path, single or all, outfolder name in you chose single
### USED FILES: 
    #  motifs_FINAL_TOTAL --> For both maize and wheat constructed previously with fields: motif_id, tf_gene_id, species, pwm_json
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/

import os
import sys
import glob
import argparse
import json
import pandas as pd
import pathlib

# GIVE ARGUMENTS TO THE SCRIPT
def parse_args(args=None):
    Description = "It changes PWMs from meme stores in json to cluster-buster"
    Epilog = "Example usage: python3 pwms_to_cluster_buster.py <PWMs_FILE_PATH>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("PWMs_FILE_PATH", help="File with all the PWMs in the json format")
    parser.add_argument(
        "-n",
        "--number",
        type=str,
        dest="number",
        default="single",
        help="Value for 'number' a file with a single or all PWMs and must be 'single' or 'all'",
        required=True
    )
    parser.add_argument(
        "-o",
        "--outfolder",
        type=str,
        dest="outfolder",
        default=None,
        help="Specify an output folder name to store the PWMs generated. Required if -n \"single\"",
    )
    
    return parser.parse_args(args)


# CONVERT A PWM FROM MEME FORMAT (stored in json format) TO CLUSTER BUSTER
def get_cluster_buster(meme_pwm_dict):
    meme_pwd_json = pd.read_json(json.dumps(meme_pwm_dict))

    cb_pwd = ''
    # lo tenemos que escalar a 100
    for index, row in meme_pwd_json.iterrows():
        a = str(int(row['A']*100))
        c = str(int(row['C']*100))
        g = str(int(row['G']*100))
        t = str(int(row['T']*100))
        new_row = '\t'.join([a, c, g, t])

        cb_pwd += new_row + '\n'
    
    return cb_pwd.rstrip()
    
    
# WRITE FILE/FILES IN CLUSTER BUSTER FORMAT PWMs
def write_cluster_buster(motif_file_path, number="single", outfolder=None):
    
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
            try:
                motif_id = fields[0]
                gene_id = fields[1]
                species = fields[2]
                pwm_dict = eval(fields[3])
            except:
                print(line)
                print(fields)

            cb_pwm = get_cluster_buster(pwm_dict)

            # IF number == single, WRITE ALL PWMs IN MOTIF FILES
            if number == "single":
                assert outfolder is not None, "Output directory name needs to be specified with -o --outfolder"

                with open(outfolder / f"{motif_id}.cb", "w") as cb:
                    cb.write(f">{motif_id}\n{cb_pwm}\n")
            
            # IF number == all, WRITE ALL PWMs IN A SINGLE MOTIF FILE
            elif number == "all" and motif_id not in motifs_in_all:
                text += f">{motif_id}\n{cb_pwm}\n"
                motifs_in_all.add(motif_id)

        if text:
            with open(FILES_DIR / f"{species}_all_motifs_cb.cb", "w") as cb:
                    cb.write(text)



def main(args=None):
    args = parse_args(args)

    write_cluster_buster(
        motif_file_path=args.PWMs_FILE_PATH,
        number=args.number,
        outfolder=args.outfolder
    )

if __name__ == '__main__':
    main()