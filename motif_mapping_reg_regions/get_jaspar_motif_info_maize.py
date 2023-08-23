### GET FILE WITH MOTIF ID, AFFECTED GENE ID, SPECIES, PWM IN JSON FORMAT FROM JASPAR DATABASE
### AND WRITE THIS INFORMATION ON THE PREVIOUS EXISTING MOTIF FILE FROM CISBP
### USED FILES: 
    # v3_to_v5_table.txt --> conversion table of Z.mays genes v3 to v5 (obtained from v3_to_v5_zmays_converter.py)
    # v4_to_v5_table.txt --> conversion table of Z.mays genes v4 to v5 (obtained from v4_to_v5_zmays_converter.py)
    # jaspar_raw.tsv --> File with motif ID, gene ID (version 3, 4 or 5) and species
    # jaspar_pwms_meme_format.txt --> File with all pwms in meme format (obtained applying jaspar2meme -bundle JASPAR2022_CORE_non-redundant_pfms_jaspar.txt)
    # motifs_FINAL_table.txt --> File with all the motif ID, gene ID, species and PWM in json format from CisBP
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/z_mays/motif_info


import pathlib
import json
from itertools import islice


PROJECT_PATH = pathlib.Path("/ngsprojects/grassgrns/data_archive/z_mays")
MOTIF_FILES_PATH = PROJECT_PATH / "motif_info"


# GET THE GENE IDs DICTIONARY WITH key = v3_genes & value = v5_genes
def get_v3_to_v5_dict(v3_to_v5_table_file):
    v3_to_v5_dict = {}
    with open(v3_to_v5_table_file, "r") as tf:
        for line in tf:
            line = line.rstrip()
            fields = line.split()
            v3_gene = fields[0]
            v5_gene = fields[1]
            v3_to_v5_dict.setdefault(v3_gene, set()).add(v5_gene) # cambiado a set

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


# GET MOTIF DICTIONARY WITH key = motif_id & value = [gene_id, species] # CAMBIAR PARA QUE TENGA EN CUENTA EL SET 
def read_all_motifs_file(motifs_file, v3_to_v5_dict, v4_to_v5_dict):
    motifs_dict = {}
    with open(motifs_file, "r") as mf:
        next(mf)
        for line in mf:
            line = line.rstrip()
            fields = line.split("\t")

            gene_name = fields[1].replace(u'\xa0', u'') # \xa0 shows up in some gene names (idk why)

            # filter if gene_id (DBID) is in the v3 genes or v4 genes and convert to v5 
            if gene_name in v3_to_v5_dict.keys():
                motif_ID = fields[0].split(".")[0]
                gene_ID_set = v3_to_v5_dict[gene_name]
                species = fields[2]
            
            elif gene_name in v4_to_v5_dict.keys():
                motif_ID = fields[0].split(".")[0]
                gene_ID_set = v4_to_v5_dict[gene_name]
                species = fields[2]
            
            else:
                motif_ID = fields[0].split(".")[0]
                gene_ID_set = {gene_name}
                species = fields[2]

            motifs_dict[motif_ID] = gene_ID_set


            # for gene_ID in gene_ID_set:
            #     motifs_dict[motif_ID] = gene_ID
                # # filter if motif_id != '.' (MOTIFS WITHOUT ID) and for species Zea mays
                # if motif_ID != "." and species == "Zea_mays":
                #     motifs_dict.setdefault(motif_ID, set()).add(tuple([gene_ID, species]))

    return motifs_dict 


# USED TO ITERATE OVER THE N LINES MAKING UP THE PWM (used in next function)
def next_n_lines(file_opened, N):
    return [x.strip() for x in islice(file_opened, N)]


# ITERATE OVER JASPAR PWMS FILE IN MEME FORMAT AND CONVERT PWM TO JSON FORMAT
def get_pwm_motif(jaspar_pwms_file, motifs_dict):
    my_lines = []
    with open(jaspar_pwms_file) as pf:
        read_next_lines = False
        for line in pf:
            line = line.rstrip()

            if line.startswith("MOTIF") and line.split()[1][:-2] in motifs_dict:
                
                motif = line.split()[1] #.split(".")[0] # esto es para quitarle el .1 al final del motif id
                gene_id_set = motifs_dict[motif.split(".")[0]]
                species =  "Zea_mays"
                read_next_lines = True
            
            elif line.startswith("letter-probability") and read_next_lines:    
                n_sites = int(line.split()[5])
                pwm = {"Pos":{}, "A":{}, "C":{}, "G":{}, "T":{}}
                index = 0
                for line in next_n_lines(pf, n_sites):
                    fields = line.split()
                    pwm["Pos"][str(index)] = index+1
                    pwm["A"][str(index)] = fields[0]
                    pwm["C"][str(index)] = fields[1]
                    pwm["G"][str(index)] = fields[2]
                    pwm["T"][str(index)] = fields[3]
                    index += 1
                
                for gene_id in gene_id_set:
                    my_lines.append((motif, gene_id, species, json.dumps(pwm, separators=(',', ':'))))
                    read_next_lines = False
    
    return my_lines
                

def main():
    v3_to_v5_table_file = MOTIF_FILES_PATH / "v3_to_v5_table.txt"
    v4_to_v5_table_file = MOTIF_FILES_PATH / "v4_to_v5_table.txt"
    jaspar_raw_file = MOTIF_FILES_PATH / "jaspar_raw.tsv"
    jaspar_pwms_file = MOTIF_FILES_PATH / "jaspar_pwms_meme_format.txt"

    v3_to_v5_dict = get_v3_to_v5_dict(v3_to_v5_table_file) 
    v4_to_v5_dict = get_v4_to_v5_dict(v4_to_v5_table_file)
    motifs_dict = read_all_motifs_file(jaspar_raw_file, v3_to_v5_dict, v4_to_v5_dict)
    
    lines = get_pwm_motif(jaspar_pwms_file, motifs_dict)

    text = ''
    for line in lines:
        text += "\t".join(line) + "\n"
    
    with open(MOTIF_FILES_PATH / "motifs_FINAL.txt", "r") as table:
        text = table.read() + "\n" + text
    
    with open(MOTIF_FILES_PATH / "motifs_FINAL_TOTAL.txt", "w") as final:
        final.write(text)
        

    # print(f"Motifs file is ready in: {MOTIF_FILES_PATH / 'motifs_FINAL_TOTAL.txt'}")
        

if __name__ == '__main__':
    main()