# Script to add ensembl IDs to the Dorothea TFTG network
#
# Input:
#   dorothea : TSV file downloaded from http://saezlab.org/tfregulons/ in May 2019 with the regulatory interactions in DoRothEA. 
#              Regulators in column labelled 'TF', and targets in column labelled 'target'. Both are gene symbols in capitals.
#   id_conversion : Human gene symbol to Ensembl ID conversion table downloaded from BioMart May 2019. Gene symbols in column 
#              "Gene.name' and Ensembl Id in column "Gene.stable.ID".
# Output: 
#   TSV file like input dorothea file but with additional columns "source_ensembl_id" and "target_ensembl_id".

library(dplyr)

# Input data
dorothea <- read.table("../input_data/prior_knowledge_networks/tfregulons_database_v01_20180216__ABC.tsv", sep = "\t", header = TRUE)
id_conversion <- read.table("../input_data/prior_knowledge_networks/mart_export_ensembl_symbol_human_100519.txt", sep = "\t", header = TRUE)

# Join the tf and target column to the BioMart table
dorothea2 <- left_join(dorothea, id_conversion, by = c("TF" = "Gene.name")) %>% dplyr::rename(source_ensembl_id = Gene.stable.ID)
dorothea3 <- left_join(dorothea2, id_conversion, by = c("target" = "Gene.name")) %>% dplyr::rename(target_ensembl_id = Gene.stable.ID)

# Remove lines where there is no conversion
dorothea4 <- dorothea3[!is.na(dorothea3$target_ensembl_id), ]
dorothea4 <- dorothea4[!is.na(dorothea4$source_ensembl_id), ]

# Extract required columns only
dorothea5 <- dorothea4 %>% select(c(source_ensembl_id, target_ensembl_id, TF, target, effect, score, is_evidence_curateddatabase, is_evidence_chipSeq, is_evidence_TFbindingMotif, is_evidence_coexpression, which_curateddatabase, which_chipSeq, which_TFbindingMotif, which_coexpression))

# Output
write.table(dorothea5, "tfregulons_database_v01_20180216__ABC_ensembl.tsv", quote = FALSE, row.names=FALSE, sep = "\t")
