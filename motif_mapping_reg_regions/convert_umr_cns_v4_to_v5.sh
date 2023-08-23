#!/bin/bash

module load liftOver/x86_64/20140423

# command used to conver bed files from v4 to v5 using a chain file obtained from (https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/assembly_chain/zea_mays/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.gz)
liftOver ./zea_mays_cnss_v4.bed /ngsprojects/grassgrns/data_archive/z_mays/genome/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain ./zea_mays_cnss_v5.bed ./unlifted/unlifted.bed
