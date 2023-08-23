### CONVERSION OF Z. MAYS GENES FROM V4 TO V5
### USED FILES: 
    # MaizeGDB_maize_pangene_2020_08.tsv --> file with all gene accessions of Z. mays in MaizeGDB
    # Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3 --> gff3 of Zea mays genome version 4
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/z_mays/motif_info

import pathlib


PROJECT_PATH = pathlib.Path("/ngsprojects/grassgrns/data_archive/z_mays")
MOTIF_FILES_PATH = PROJECT_PATH / "motif_info"

# GET A SET OF V4 GENES IN Zea mays (NON-USED)
# def get_v4_geneset(gff_file):
#     v4_genes = set()

#     with open(gff_file, "r") as gff:
#         for line in gff:
#             line = line.rstrip()
#             if not line.startswith('#'):
#                 fields = line.split()
#                 line_type = fields[2]

#                 if line_type == "gene":
#                     gene_ID = fields[8].split(";")[0].split(":")[1]
#                     if gene_ID not in v4_genes:
#                         v4_genes.add(gene_ID)

#     return v4_genes


# GET DICTIONARY WHERE V4 GENES ARE KEY AND LIST OF V5 GENES ARE THE VALUE
# def get_v4_to_v5_dict(tsv_file, v4_genes): (NON-USED)
def get_v4_to_v5_dict(tsv_file):
    v4_to_v5_dict = {}
    total_v5 = []
    
    with open(tsv_file, "r") as tsv:
        for line in tsv:
            v4_set = set()
            v5_set = set()
            
            line = line.rstrip()
            fields = line.split()

            for field in fields:
                if 'Zm00001d' in field:
                # if field in v4_genes: (NON-USED)
                    v4_set.add(field)

                elif field.startswith('Zm00001eb'):
                    v5_set.add(field)
                    total_v5.append(field)

            if v4_set and v5_set:
                for gene_4 in v4_set:
                    for gene_5 in v5_set:
                        v4_to_v5_dict.setdefault(gene_4, set()).add(gene_5)


    return v4_to_v5_dict
    
    

def main():
    # v4_gff_file = MOTIF_FILES_PATH / "Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3" (NON-USED)
    tsv_file = MOTIF_FILES_PATH / "MaizeGDB_B73_pangene_2020_11.tsv"

    # v4_genes = get_v4_geneset(v4_gff_file) (NON-USED)
    # v4_to_v5_dict = get_v4_to_v5_dict(tsv_file, v4_genes) (NON-USED)
    v4_to_v5_dict = get_v4_to_v5_dict(tsv_file)

    # WRITE A FILE WITH 1st COLUMN BEING V4 GENES AND 2nd COLUMN BEING A LIST OF V5 GENES
    with open(MOTIF_FILES_PATH / "v4_to_v5_table.txt", "w") as test:
        text = ""
        for key, values in v4_to_v5_dict.items():
            if values:
                for value in values:
                    text += key + "\t" + value + "\n"
        
        test.write(text)
    
    print(f"v4 to v5 file is ready in: {MOTIF_FILES_PATH / 'v4_to_v5_table.txt'}")
    

if __name__ == '__main__':
    main()