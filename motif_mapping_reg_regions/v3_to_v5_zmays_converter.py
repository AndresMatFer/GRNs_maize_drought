### CONVERSION OF Z. MAYS GENES FROM V3 TO V5
### USED FILES: 
    # MaizeGDB_maize_pangene_2020_08.tsv --> file with all gene accessions of Z. mays in MaizeGDB
    # Zea_mays.AGPv3.22.gff3 --> gff3 of Zea mays genome version 3
### ALL FILES CAN BE FOUND IN: /ngsprojects/grassgrns/data_archive/z_mays/motif_info

import pathlib


PROJECT_PATH = pathlib.Path("/ngsprojects/grassgrns/data_archive/z_mays")
MOTIF_FILES_PATH = PROJECT_PATH / "motif_info"

# GET A SET OF V3 GENES IN Zea mays
def get_v3_geneset(gff_file):
    v3_genes = set()

    with open(gff_file, "r") as gff:
        for line in gff:
            line = line.rstrip()
            if not line.startswith('#'):
                fields = line.split()
                line_type = fields[2]

                if line_type == "gene":
                    gene_ID = fields[8].split(";")[0].split(":")[1]
                    if gene_ID not in v3_genes:
                        v3_genes.add(gene_ID)

    return v3_genes


# GET DICTIONARY WHERE V3 GENES ARE KEY AND LIST OF V5 GENES ARE THE VALUE
def get_v3_to_v5_dict(tsv_file, v3_genes):
    v3_to_v5_dict = {}
    total_v5 = []
    
    with open(tsv_file, "r") as tsv:
        for line in tsv:
            v3_set = set()
            v5_set = set()
            
            line = line.rstrip()
            fields = line.split()

            for field in fields:
                # there are genes v3 like (AC187090.3_FG011)
                if field.startswith("B73v3") and "_".join(field.split("_")[1:]) in v3_genes:
                    v3_gene_name = "_".join(field.split("_")[1:])
                    v3_set.add(v3_gene_name)

                elif field.startswith('Zm00001eb'):
                    v5_gene_name = field
                    v5_set.add(field)
                    total_v5.append(field)
            
            # if sets are empty (because there are no v3 or v5)
            if v3_set and v5_set:
                for gene_3 in v3_set:
                    for gene_5 in v5_set:
                        # v3_to_v5_dict[gene] = v5_list
                        v3_to_v5_dict.setdefault(gene_3, set()).add(gene_5)


    return v3_to_v5_dict
    
    

def main():
    v3_gff_file = MOTIF_FILES_PATH / "Zea_mays.AGPv3.22.gff3"
    tsv_file = MOTIF_FILES_PATH / "MaizeGDB_B73_pangene_2020_11.tsv"

    v3_genes = get_v3_geneset(v3_gff_file)
    v3_to_v5_dict = get_v3_to_v5_dict(tsv_file, v3_genes)

    # WRITE A FILE WITH 1st COLUMN BEING V3 GENES AND 2nd COLUMN BEING A LIST OF V5 GENES
    with open(MOTIF_FILES_PATH / "v3_to_v5_table.txt", "w") as test:
        text = ""
        for key, values in v3_to_v5_dict.items():
            if values:
                for value in values:
                    text += key + "\t" + value + "\n"
        
        test.write(text)
    
    print(f"v3 to v5 file is ready in: {MOTIF_FILES_PATH / 'v3_to_v5_table.txt'}")
    

if __name__ == '__main__':
    main()