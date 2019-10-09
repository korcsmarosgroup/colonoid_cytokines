# Functional analysis of DEGs found in the biopsy data but not in the colonoid data (pre-network generation)
#
# Input: 
#       all_org_files - File paths for the preprocessed (during 'network_generation.py') differential expression data for each 
#                       of the cytokine treated organoids. With adj p value filter applied.
#       biop_file - Excel file of preprocessed biopsy differential expression (provided by Polchronis)
#
# Output:
#     Text file for each type of functional enrichment with all significant results for all combinations of biopsy uc/cd/cytokine overlaps: Reactome, KEGG, GO.
#     Venn diagram of the DEGs shared (or not) between the biopsies UC/CD and the colonoids.

##### Set up #####

library(tidyr)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(dplyr) # Must load after Biomart
library(VennDiagram)

setwd("../input_data/differential_expression_datasets/")
all_org_files <- list.files(pattern = "_degs_adjp_0.01.txt$")
biop_file <- "Biopsies/GSE16879DEGs/GSE16879DEGsUCandcCDvsHCGEO2R.xlsx"

lfc <- TRUE # true or false to apply the lfc filter 1.5

##### Process organoid DEGs #####

# Import deg lists
idx <- 0
for (file in all_org_files) {
  
  idx <- idx+1
  
  # Open file
  degs <- read.delim(file, sep = "\t",header=TRUE,na.strings=c("","NA"))
  
  if (lfc) {
    # OR INSTEAD Remove columns, filter for lfc >= 1.5
    degs <- degs  %>% filter((logFC >1.5) | (logFC < -1.5)) %>% dplyr::select(ENSEMBL.ID, SYMBOL)
  } else {
    # Remove columns
    degs <- degs %>% dplyr::select(ENSEMBL.ID, SYMBOL)
    
  }
    
  # Get cytokine name
  cytokine <- str_replace(file, "_degs_adjp_0.01.txt", "")
  degs <- degs %>% mutate(cytokine = cytokine)
  
  #join degs together
  if (idx == 1){
    all_degs <- degs
  } else {
    all_degs <- full_join(all_degs, degs, by=c("ENSEMBL.ID", "SYMBOL"))
  }
  
}

# Join the cytokine columns together
all_degs2 <- unite(all_degs, cytokines, 3:6, sep="; ") %>%
  mutate(cytokines = str_replace_all(cytokines, '; NA,?','')) %>%
  mutate(cytokines = str_replace_all(cytokines, 'NA,?; ',''))

##### Process biopsy DEGs #####

#GEO file
cd_degs <- read_excel(biop_file, sheet ="cCDvsHC")
uc_degs <- read_excel(biop_file, sheet ="UCvsHC")
# Seperate gene symbols
cd_degs <- separate_rows(cd_degs,Gene.symbol,sep="///")
uc_degs <- separate_rows(uc_degs,Gene.symbol,sep="///")

if (lfc) {
  # OR INSTEAD Filter for p adj and lfc 1.5
  cd_degs <- cd_degs %>% filter(adj.P.Val <= 0.01 & (logFC >= 1.5 | logFC <= -1.5)) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_cd = "TRUE") %>% unique() %>% drop_na()
  uc_degs <- uc_degs %>% filter(adj.P.Val <= 0.01 & (logFC >= 1.5 | logFC <= -1.5)) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_uc = "TRUE") %>% unique() %>% drop_na()
} else {
  # Filter for p adj <= 0.01, get gene names
  cd_degs <- cd_degs %>% filter(adj.P.Val <= 0.01) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_cd = "TRUE") %>% unique() %>% drop_na()
  uc_degs <- uc_degs %>% filter(adj.P.Val <= 0.01) %>% dplyr::select(c(Gene.symbol)) %>% mutate(in_uc = "TRUE") %>% unique() %>% drop_na()
}
  
# Convert gene symbol to ensembl ID
row.names(cd_degs) <- cd_degs$Gene.symbol
cd_degs$Ensembl <- mapIds(org.Hs.eg.db, keys = row.names(cd_degs), keytype = "SYMBOL", column="ENSEMBL",multiVals="first")
row.names(uc_degs) <- uc_degs$Gene.symbol
uc_degs$Ensembl <- mapIds(org.Hs.eg.db, keys = row.names(uc_degs), keytype = "SYMBOL", column="ENSEMBL",multiVals="first")
cd_degs <- cd_degs %>% na.omit() %>% dplyr::select(c(Ensembl, in_cd)) %>% unique()
uc_degs <- uc_degs %>% na.omit() %>% dplyr::select(c(Ensembl, in_uc)) %>% unique()

##### Extract biopsy DEGs not in colonoid data ######

# DEGs in all 3 datasets
cd_uc_colonoids <- uc_degs %>% filter((Ensembl %in% cd_degs$Ensembl) & (Ensembl %in% all_degs2$ENSEMBL.ID)) %>% mutate(list = "cd_and_uc_and_colonoid")

# DEGs shared between uc and cd datasets
uc_cd <- uc_degs %>% filter(Ensembl %in% cd_degs$Ensembl) %>% mutate(list = "cd_and_uc_all")
uc_cd_notcolonoid <- uc_degs %>% filter((Ensembl %in% cd_degs$Ensembl) & (!Ensembl %in% all_degs2$ENSEMBL.ID)) %>% mutate(list = "cd_and_uc_notcolonoid")

# DEGs in colonoids and UC and visa versa
uc_colonoid_shared <- uc_degs %>% filter(Ensembl %in% all_degs2$ENSEMBL.ID) %>% mutate(list = "uc_and_colonoids_all")
cd_colonoid_shared <- cd_degs %>% filter(Ensembl %in% all_degs2$ENSEMBL.ID) %>% mutate(list = "cd_and_colonoids_all")
uc_colonoid_notcd <- uc_degs %>% filter((Ensembl %in% all_degs2$ENSEMBL.ID)&(!Ensembl %in% cd_degs$Ensembl)) %>% mutate(list = "uc_and_colonoids_notcd")
cd_colonoid_notuc <- cd_degs %>% filter((Ensembl %in% all_degs2$ENSEMBL.ID)&(!Ensembl %in% uc_degs$Ensembl)) %>% mutate(list = "cd_and_colonoids_notuc")

cd_only <- cd_degs %>% filter(!Ensembl %in% all_degs2$ENSEMBL.ID) %>% mutate(list = "cd_notcolonoid_all")
uc_degs_filt <- uc_degs %>% filter(!Ensembl %in% all_degs2$ENSEMBL.ID)  %>% mutate(list = "uc_notcolonoid_all")

cd_spec_filt <- cd_degs %>% filter((!Ensembl %in% uc_degs$Ensembl)&(!Ensembl %in% all_degs2$ENSEMBL.ID))  %>% mutate(list = "cd_notcolonoid_notuc")
uc_spec_filt <- uc_degs %>% filter((!Ensembl %in% cd_degs$Ensembl)&(!Ensembl %in% all_degs2$ENSEMBL.ID))   %>% mutate(list = "uc_notcolonoid_notcd")

colonoids_only <- all_degs2 %>% filter((!ENSEMBL.ID %in% cd_degs$Ensembl)&(!ENSEMBL.ID %in% uc_degs$Ensembl)) %>% mutate(list = "colonoids_notuc_notcd") %>% dplyr::rename(Ensembl = ENSEMBL.ID)
  
input_dfs <- list(cd_uc_colonoids, uc_cd, uc_cd_notcolonoid, uc_colonoid_shared,cd_colonoid_shared,uc_colonoid_notcd,cd_colonoid_notuc,cd_only, uc_degs_filt,cd_spec_filt,uc_spec_filt,colonoids_only)

##### Carry out functional analysis ######

# Empty df's for data
all_reactome <- data.frame()
all_kegg <- data.frame()
all_go <- data.frame()

for (input in input_dfs) {
  
  #convert to entrez ids
  input_ent <- bitr(input$Ensembl, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #Get list name
  input_name <- input$list[[1]]
  
  #### REACTOME ####
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=input_ent$ENTREZID,pvalueCutoff=0.05,readable=T)
  #Get reactome results
  genes_reactome2 <- as.data.frame(genes_reactome)
  genes_reactome2 <- mutate(genes_reactome2, group = input_name)
  
  #Append reactome results to previous results
  all_reactome <- rbind(all_reactome, genes_reactome2)
  
  #### KEGG ####
  
  #KEGG pathway enrichment
  genes_kegg <- enrichKEGG(gene=input_ent$ENTREZID ,organism = 'hsa',pvalueCutoff = 0.05)
  #Get KEGG results
  genes_kegg2 <- as.data.frame(genes_kegg)
  genes_kegg2 <- mutate(genes_kegg2, group = input_name)
  
  #Append reactome results to previous results
  all_kegg <- rbind(all_kegg, genes_kegg2)
  
  #### GO ####
  
  #GO over representation
  genes_go <- enrichGO(gene=input_ent$ENTREZID,OrgDb = "org.Hs.eg.db", ont= "BP", pvalueCutoff=0.05)
  #Get GO results
  genes_go2 <- as.data.frame(genes_go)
  genes_go2 <- mutate(genes_go2, group = input_name)
  
  #Append reactome results to previous results
  all_go <- rbind(all_go, genes_go2)
}
  
# Save outputs

write.table(all_reactome, "../biopsy_colonoid_functionalEnrichment_Reactome_0.05_geo_lfc1.5.txt",quote = F, sep = "\t", row.names = F)
write.table(all_kegg, "../biopsy_colonoid_functionalEnrichment_KEGG_0.05_geo_lfc1.5.txt",quote = F, sep = "\t", row.names = F)
write.table(all_go, "../biopsy_colonoid_functionalEnrichment_GO_0.05_geo_lfc1.5.txt",quote = F, sep = "\t", row.names = F)

##### Venn diagram ####

pdf("../biopsy_colonoid_shared_DEGs_EnsemblIDs_qval0.01_lfc1.5.pdf",width=7,height=7) 
grid.newpage()
draw.triple.venn(area1 = nrow(all_degs2), area2 = nrow(uc_degs), area3 = nrow(cd_degs), n12 = nrow(uc_colonoid_shared), n23 = nrow(uc_cd), n13 = nrow(cd_colonoid_shared), 
                 n123 = nrow(cd_uc_colonoids), category = c("Colonoid DEGs", "UC DEGs", "CD DEGs"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()
