# Script to generate network to visualise the cyto-tf-deg interaction for the colonoid-cytokine data
#
# Author: Agatha Treveil
# Date: May 2020
# 
# Input: Dataframe of all tf-deg interactions with cytokine categories - generated previously
#

##### Set up #####
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(stringr)

setwd("/Users/treveila/Documents/cytokines_colonoids/analysis/tf-deg_investigation/tf-deg_ints/")
whole_network <- read.csv("tf-deg_all_ints.txt", sep = "\t")

outfile <- "il13_cyto-tf-deg_summary_network.txt"

cyto <- "il13"

##### Cyto-TF network #####

# Map to labels we want
whole_network$deg_category <- plyr::mapvalues(whole_network$deg_category, from=c("ifng,il13,tnfa", "ifng,il13","ifng,il13,il17,tnfa","il17,tnfa","ifng,il13,il17", "ifng,il17,tnfa"),
                                    to=c("il13_ifng_tnfa", "il13_ifng","il13_ifng_tnfa_il17","tnfa_il17","il13_ifng_il17","ifng_tnfa_il17")) %>% str_replace_all(",","_")
whole_network$tf_category <- plyr::mapvalues(whole_network$tf_category, from=c("ifng,il13,tnfa", "ifng,il13","ifng,il13,il17,tnfa","il17,tnfa","ifng,il13,il17","ifng,il17,tnfa"),
                                   to=c("il13_ifng_tnfa", "il13_ifng","il13_ifng_tnfa_il17","tnfa_il17","il13_ifng_il17","ifng_tnfa_il17")) %>% str_replace_all(",","_")

# Filter tf_categries containing cyto and count (cytokine specific)
tfs_sp <- whole_network %>% select(tf_ensembl_id, tf_category) %>% unique() %>% 
  filter(tf_category == cyto) %>% group_by(tf_category) %>% tally()

# Transform
rownames(tfs_sp) <- tfs_sp$tf_category
tfs_sp <- tfs_sp %>% select(n) %>% t() 

# Filter tf_categries containing cyto and count (not cytokine specific)
tfs_sh <- whole_network %>% select(tf_ensembl_id, tf_category) %>% unique() %>% 
  filter(grepl(cyto, tf_category) & (tf_category != cyto)) %>% group_by(tf_category) %>% tally()

# Transform
rownames(tfs_sh) <- tfs_sh$tf_category
tfs_sh <- tfs_sh %>% select(n) %>% t() 

# Create rows
edges <- data.frame()
row <- data.frame(source = paste0("cyto_",cyto), target = paste0("tf_",cyto,"_sp"), target_size = sum(tfs_sp[1,]), stringsAsFactors=FALSE)
row <- cbind(row,tfs_sp)
edges <- bind_rows(edges,row)  

row <- data.frame(source = paste0("cyto_",cyto), target = paste0("tf_",cyto,"_sh"), target_size = sum(tfs_sh[1,]), stringsAsFactors=FALSE)
row <- cbind(row,tfs_sh)
edges <- bind_rows(edges,row)  

##### Cytokine TFs #####

# Get cytokine DEGs (specific and shared)
deg <- whole_network %>% filter(grepl(cyto, deg_category)) %>% select(deg_ensembl_id, deg_category, tf_category) %>% unique()

# Get the tf categories for these degs
#deg_tf <- deg %>% select(deg_ensembl_id,tf_category) %>% unique()
deg_tf <- spread(deg,tf_category, tf_category, convert=T)
deg_tf <- unite(deg_tf, tf_cat, 3:ncol(deg_tf), sep = ";", na.rm = T)


##### Get DEGS with only cytokine specific tfs #####

# Get cytokine targeted degs where tf is cytokine specific
deg_tf_sp <- deg_tf %>% filter(tf_cat == cyto) %>% unique()

# Get deg categories for these
deg_tf_sp2 <- deg_tf_sp %>% select(deg_category) %>% group_by(deg_category) %>% tally()
rownames(deg_tf_sp2) <- deg_tf_sp2$deg_category
deg_tf_sp2 <- deg_tf_sp2 %>% select(n) %>% t() 

# Add row
row <- data.frame(source = paste0("tf_",cyto,"_sp"), target=paste0("deg_", cyto, "_sp"), target_size =nrow(deg_tf_sp))
row <-cbind(row,deg_tf_sp2)
edges <- bind_rows(edges,row)  

##### Get DEGs which have both cytokine specific and shared TFs #####

# Get cytokine targeted degs where tfs are cytokine specific and not
deg_tf_spsh <- deg_tf %>% filter(grepl(paste0(";",cyto,";"), tf_cat) | grepl(paste0("^",cyto,";"), tf_cat) | grepl(paste0(";",cyto,"$"), tf_cat)) %>% unique()

# Get deg categories for these
deg_tf_spsh2 <- deg_tf_spsh %>% select(deg_category) %>% group_by(deg_category) %>% tally()
rownames(deg_tf_spsh2) <- deg_tf_spsh2$deg_category
deg_tf_spsh2 <- deg_tf_spsh2 %>% select(n) %>% t() 

# Add rows
row <- data.frame(source = paste0("tf_",cyto,"_sp"), target=paste0("deg_", cyto, "_spsh"), target_size =nrow(deg_tf_spsh))
row <-cbind(row,deg_tf_spsh2)
edges <- bind_rows(edges,row)  
row <- data.frame(source = paste0("tf_",cyto,"_sh"), target=paste0("deg_", cyto, "_spsh"), target_size =nrow(deg_tf_spsh))
row <-cbind(row,deg_tf_spsh2)
edges <- bind_rows(edges,row)  

##### Get DEGs which have only cytokine shared TFs #####

# Get cytokine targeted degs where tfs are all cytokine shared
deg_tf_sh <- deg_tf %>% filter(!grepl(paste0(";",cyto,";"), tf_cat) & !grepl(paste0("^",cyto,";"), tf_cat) & !grepl(paste0(";",cyto,"$"), tf_cat) & (tf_cat != cyto)) %>% unique()

# Get deg categories for these
deg_tf_sh2 <- deg_tf_sh %>% select(deg_category) %>% group_by(deg_category) %>% tally()
rownames(deg_tf_sh2) <- deg_tf_sh2$deg_category
deg_tf_sh2 <- deg_tf_sh2 %>% select(n) %>% t() 

# Add rows
row <- data.frame(source = paste0("tf_",cyto,"_sh"), target=paste0("deg_", cyto, "_sh"), target_size =nrow(deg_tf_sh))
row <-cbind(row,deg_tf_sh2)
edges <- bind_rows(edges,row) 

##### Save network file ####

# Convert all NA to 0
edges[is.na(edges)] <- 0

write.table(edges, file = outfile, quote = F, row.names = F, sep = "\t")
