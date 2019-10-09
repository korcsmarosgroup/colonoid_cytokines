# Functional analysis on the DEGs targeted by each cytokine (grouping by which combination 
# of cytokine they are targeted by based on the causal networks.
#
# Input: input - Path to cytokine-deg direct network generated in the script Cytokine_to_tf_or_deg_direct.py.
#
# Output: Text file for each type of functional enrichment with all significant results for all subnetworks: Reactome, KEGG, GO.

##### Set up #####

# Packages
library(dplyr)
library(ReactomePA)
library(clusterProfiler)
library("org.Hs.eg.db")

# Input files
input <- read.delim("../input_data/causal_networks/cytokine_deg_networks/cytokine_deg_direct/AllCytokines_to_deg_direct_network.txt", header=T)

# Loop through groups of gene id based on cytokine_all column
input_split <- split(input, f = input$cytokine_all)

# Empty df's for data
all_reactome <- data.frame()
all_kegg <- data.frame()
all_go <- data.frame()

##### Functional enrichment analysis #####

for (list in input_split) {
  
  #convert to entrez ids
  input_e <- bitr(list$deg_ensembl_id, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #Get list name
  list_f <- list$cytokine_all[1]
  
  #### REACTOME ####
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=input_e$ENTREZID,pvalueCutoff=0.1,readable=T)
  #Get reactome results
  genes_reactome2 <- as.data.frame(genes_reactome)
  genes_reactome2 <- mutate(genes_reactome2, group = list_f)
  
  #Append reactome results to previous results
  all_reactome <- rbind(all_reactome, genes_reactome2)
  
  #### KEGG ####
  
  #KEGG pathway enrichment
  genes_kegg <- enrichKEGG(gene=input_e$ENTREZID ,organism = 'hsa',pvalueCutoff = 0.05)
  #Get KEGG results
  genes_kegg2 <- as.data.frame(genes_kegg)
  genes_kegg2 <- mutate(genes_kegg2, group = list_f)
  
  #Append reactome results to previous results
  all_kegg <- rbind(all_kegg, genes_kegg2)
  
  #### GO ####
  
  #GO over representation
  genes_go <- enrichGO(gene=input_e$ENTREZID,OrgDb = "org.Hs.eg.db", ont= "BP", pvalueCutoff=0.05)
  #Get GO results
  genes_go2 <- as.data.frame(genes_go)
  genes_go2 <- mutate(genes_go2, group = list_f)
  
  #Append reactome results to previous results
  all_go <- rbind(all_go, genes_go2)
}

##### Save outputs #####

write.table(all_reactome, "FunctionalEnrichment_Reactome_0.1.txt",quote = F, sep = "\t", row.names = F)
write.table(all_kegg, "FunctionalEnrichment_KEGG_0.05.txt",quote = F, sep = "\t", row.names = F)
write.table(all_go, "FunctionalEnrichment_GO_0.05.txt",quote = F, sep = "\t", row.names = F)
