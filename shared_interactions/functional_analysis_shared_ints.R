# Functional analysis on shared interacctions (input list of human ensembl ids)
#
# Input: network - All interactions which are shared between cytokine networks with the 'Cytokines' column indicating which networks they are present in. 
#                 Created in the combine_network_data.py script (and edited to remove interactions in only one network).
#
# Output: Text file for each type of functional enrichment with all significant results for all subnetworks: Reactome, KEGG, GO.

##### Set up #####

# Packages
library(dplyr)
library(ReactomePA)
library(clusterProfiler)

#Input files
input <- read.delim("../input_data/causal_networks/cytokine_deg_networks/shared_interactions/Interactions_2ormore.txt", header=T)

#Loop through groups of gene id based on cytokine_all column
input_split <- split(input, f = input$Cytokines)

# Empty df's for data
all_reactome <- data.frame()
all_kegg <- data.frame()
all_go <- data.frame()

##### Functional analysis #####

for (list in input_split) {
  
  # get all ensembl ids in one list
  id_list <- unlist(list[c("source", "target")], use.names = FALSE)
  id_list <- id_list[!duplicated(id_list)]
  
  #convert to entrez ids
  input_e <- bitr(id_list, fromType='ENSEMBL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #Get list name
  list_f <- list$Cytokines[1]
  
  #### REACTOME ####
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=input_e$ENTREZID,qvalueCutoff=0.1,readable=T)
  #Get reactome results
  genes_reactome2 <- as.data.frame(genes_reactome)
  genes_reactome2 <- mutate(genes_reactome2, group = list_f)
  
  #Append reactome results to previous results
  all_reactome <- rbind(all_reactome, genes_reactome2)
  
  #### KEGG ####
  
  #KEGG pathway enrichment
  genes_kegg <- enrichKEGG(gene=input_e$ENTREZID ,organism = 'hsa',qvalueCutoff = 0.1)
  #Get KEGG results
  genes_kegg2 <- as.data.frame(genes_kegg)
  genes_kegg2 <- mutate(genes_kegg2, group = list_f)
  
  #Append reactome results to previous results
  all_kegg <- rbind(all_kegg, genes_kegg2)
  
  #### GO ####
  
  #GO over representation
  genes_go <- enrichGO(gene=input_e$ENTREZID,OrgDb = "org.Hs.eg.db", ont= "BP", qvalueCutoff=0.1)
  # Remove redundancy of GO terms
  genes_go1 <- simplify(genes_go, cutoff=0.7, select_fun=min)
  #Get GO results
  genes_go2 <- as.data.frame(genes_go1)
  genes_go2 <- mutate(genes_go2, group = list_f)
  
  #Append reactome results to previous results
  all_go <- rbind(all_go, genes_go2)
}

##### Save outputs #####

write.table(all_reactome, "FunctionalEnrichment_Reactome_0.1.txt",quote = F, sep = "\t", row.names = F)
write.table(all_kegg, "FunctionalEnrichment_KEGG_0.1.txt",quote = F, sep = "\t", row.names = F)
write.table(all_go, "FunctionalEnrichment_GO_0.1.txt",quote = F, sep = "\t", row.names = F)
