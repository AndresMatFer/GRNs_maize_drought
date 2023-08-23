### GET FILE WITH MOTIF ID, AFFECTED GENE ID, SPECIES, PWM IN JSON FORMAT
### USED FILES:
    # v3_to_v5_table.txt --> conversion table of Z.mays genes v3 to v5 (obtained from v3_to_v5_zmays_converter.py)
    # TF_Information_all_motifs.txt --> File with motif data for all TFs in Z. mays
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/z_mays/motif_info

import pathlib
import os
import pandas as pd


PROJECT_PATH = pathlib.Path("/ngsprojects/grassgrns/data_archive/z_mays/")
MOTIF_FILES_PATH = PROJECT_PATH / "motif_info"
PWM_FILES = MOTIF_FILES_PATH / "pwms_all_motifs"

NICO_TEST_FILE = "/scratch/tmp/niman/maizeGDB_GRAMENE_ensemble_geneid_mappings_onlygffgenes.txt"

# GET THE GENE IDs DICTIONARY WITH key = v3_genes & value = v5_genes
def get_v3_to_v5_dict(v3_to_v5_table_file):
    v3_to_v5_dict = {}
    with open(v3_to_v5_table_file, "r") as tf:
        for line in tf:
            line = line.rstrip()
            fields = line.split()
            v3_gene = fields[0]
            v5_gene = fields[1]
            v3_to_v5_dict.setdefault(v3_gene, set()).add(v5_gene) # +cambiar a set de v5 (puede haver varios v5 asociados a un v3)

    return v3_to_v5_dict


# GET THE GENE IDs DICTIONARY WITH key = v4_genes & value = v5_genes
def get_v4_to_v5_dict(v4_to_v5_table_file):
    v4_to_v5_dict = {}
    with open(v4_to_v5_table_file, "r") as tf:
        for line in tf:
            line = line.rstrip()
            fields = line.split()
            v4_gene = fields[0]
            v5_gene = fields[1]
            v4_to_v5_dict.setdefault(v4_gene, set()).add(v5_gene) # cambiado a set

    return v4_to_v5_dict


# GET cisBP uniprot dictionary
def get_cisbp_to_plaza_dict(cisbp_plaza_file, v4_to_v5_dict):
    cisbp_plaza_dict = {}
    with open(cisbp_plaza_file) as cp:
        for line in cp:
            line = line.rstrip()
            fields = line.split()
            cisbp = fields[0]
            plaza = fields[1]
            if plaza in v4_to_v5_dict:
                cisbp_plaza_dict[cisbp] = v4_to_v5_dict[plaza].pop()

    return cisbp_plaza_dict

# GET MOTIF DICTIONARY WITH key = motif_id & value = [gene_id, species]
def read_all_motifs_file(motifs_file, v3_to_v5_dict, cisbp_plaza_dict):
    motifs_dict = {}
    with open(motifs_file, "r") as mf:
        next(mf)
        for line in mf:
            line = line.rstrip()
            fields = line.split()

            # filter if gene_id (DBID) is in the v3 genes and convert to v5
            if fields[5] in v3_to_v5_dict:
                motif_ID = fields[3]
                gene_ID_set = v3_to_v5_dict[fields[5]]
                species = fields[7]

                # filter if motif_id != '.' (MOTIFS WITHOUT ID) and for species Zea mays
                if motif_ID != "." and species == "Zea_mays":
                    for gene_ID in gene_ID_set:
                        motifs_dict.setdefault(motif_ID, set()).add(tuple([gene_ID, species]))

            elif "_MAIZE" in fields[5]:
                motif_ID = fields[3]
                try:
                    gene_ID = cisbp_plaza_dict[fields[5]]
                except KeyError: # there are two genes from a different line (we can ignore them)
                    print(fields[5]) # just print them so we can be sure that we are only missing these two
                species = fields[7]

                # filter if motif_id != '.' (MOTIFS WITHOUT ID) and for species Zea mays
                if motif_ID != "." and species == "Zea_mays":
                    if motif_ID not in motifs_dict:
                        motifs_dict[motif_ID] = set()
                    motifs_dict[motif_ID].add(tuple([gene_ID, species]))


    return motifs_dict


# return the PWM of a motif as a json or return false if the PWM is empty
def get_pwm_info(pwm_path_file):
    if os.stat(pwm_path_file).st_size != 12: # 12 is the size for empty pwm files, with only: (Pos	A	C	G	T   )
        pwm_datafame = pd.read_csv(pwm_path_file, sep='\t')
        pwm_json = pwm_datafame.to_json()

        return pwm_json

    return False


def main():
    # in CisBP TFs are in v3 version, we need to update to v5 with the conversion file
    v3_to_v5_table_file = MOTIF_FILES_PATH / "v3_to_v5_table.txt"
    v4_to_v5_table_file = MOTIF_FILES_PATH / "v4_to_v5_table.txt"
    motifs_file = MOTIF_FILES_PATH / "TF_Information_all_motifs.txt" # not use "TF_Information_all_motifs_plus.txt"
    cisbp_plaza_file = MOTIF_FILES_PATH / "Uniprot_PLAZA_GeneID_Mapping_noPH207.txt"

    v3_to_v5_dict = get_v3_to_v5_dict(v3_to_v5_table_file) ## <-- we can check here with NICO_TEST_FILE
    v4_to_v5_dict = get_v4_to_v5_dict(v4_to_v5_table_file)
    cisbp_plaza_dict = get_cisbp_to_plaza_dict(cisbp_plaza_file, v4_to_v5_dict)
    motifs_dict = read_all_motifs_file(motifs_file, v3_to_v5_dict, cisbp_plaza_dict)

    text = ''

    for motif, genes in motifs_dict.items():
        pwm_path = PWM_FILES / f"{motif}.txt"
        pwm_json = get_pwm_info(pwm_path)

        if pwm_json:
            for gene_and_species in genes:
                # gene_and_species = genes.pop()
                text += '\t'.join([motif.split('_')[0]] + list(gene_and_species) + [pwm_json]) + '\n'

    # write the result table into a file
    with open(MOTIF_FILES_PATH / "motifs_FINAL.txt", "w") as table:
        table.write("motif_id\tgene_id\tspecies\tPWM\n" + text)

    print(f"Motifs file is ready in: {MOTIF_FILES_PATH / 'motifs_FINAL.txt'}")


if __name__ == '__main__':
    main()