# Script to visualise the cytokine networks as chord diagrams: cytokine-tf, cytokine-deg and tf-deg
#
# Input: 
#     working_dir - Must contain all cytokine-deg networks to be visualised. Network file names must end with 'omnip_deg_network.txt'.
#     symbol_convert - Biomart downloaded conversion file from human Ensembl IDs to gene symbols.
#     colours_file - Text file containing the colour (hex code) for each cytokine group eg. 'Il13,tnfa'. Group in column 'category' and hex code in column 'colour'.
# Output:
#     The script creates 3 chord diagrams (cyto-tf, cyto-deg and tf-deg) but it doesn't save them automatically.

#### SETUP ####

library(circlize)
library(dplyr)
library(stringr)

setwd("../input_data/causal_networks/cytokine_deg_networks/")
symbol_convert = "../input_data/id_conversion/mart_export_ensembl_symbol_human_100519.txt" 
colours_file = "../input_data/plot_colours_all.txt"

# Import data
all_files <- list.files(pattern = "omnip_deg_network.txt$")
colours = read.delim(colours_file)

#### GET ALL NETWORK DATA ####

idx <- 0
for (file in all_files) {
  
  idx <- idx+1
  
  # Open file
  network <- read.delim(file, sep = "\t",header=TRUE,na.strings=c("","NA"))
  
  # Get cytokine name
  cytokine <- str_replace(file, "_omnip_deg_network.txt", "")
    
  # Remove columns
  network2 <- network %>% select(level, source_ensembl_id, target_ensembl_id, source_symbol, target_symbol)
  
  # Add columns - level and category
  network2 <- network2 %>% mutate(source_layer = if_else(level == "tf-target", "tf", 
                                                        if_else(grepl('receptor-interm1',level), "receptor",
                                                                if_else(grepl('interm-interm',level), "signalling",
                                                                        if_else(grepl('final_interm-tf',level), "signalling",
                                                                                "cytokine")))))
  network2 <- network2 %>% mutate(target_layer = if_else(level == "tf-target", "target", 
                                                         if_else(grepl('receptor-interm1',level), "signalling",
                                                                 if_else(grepl('interm-interm',level), "signalling",
                                                                         if_else(grepl('final_interm-tf',level), "tf",
                                                                                 "receptor")))))
  network2 <- network2 %>% mutate(source_category = cytokine)
  network2 <- network2 %>% mutate(target_category = cytokine)
  
  #Get all tfs
  tfs <- network2 %>% select(tf_ensembl_id = source_ensembl_id, tf_symbol = source_symbol, layer = level, cytokine = source_category) # target_ensembl_id, target_symbol
  tfs <- tfs %>% filter(str_detect(layer, "tf-target")) %>% select(-layer) %>% unique() #(layer, "final_interm-tf")
  
  #join dfs together
  if (idx == 1){
    all_tfs <- tfs
  } else {
    all_tfs <- rbind(all_tfs, tfs)
  }
  
  #Get all degs
  degs <- network2 %>% select(deg_ensembl_id = target_ensembl_id, deg_symbol = target_symbol, layer = level, cytokine = source_category)
  degs <- degs %>% filter(str_detect(layer, "tf-target")) %>% select(-layer) %>% unique()
  
  #join degs together
  if (idx == 1){
    all_degs <- degs
  } else {
    all_degs <- rbind(all_degs, degs)
  }
  
  #Get all tf-degs
  tf_deg <- network2 %>% select(deg_ensembl_id = target_ensembl_id, deg_symbol = target_symbol, tf_ensembl_id = source_ensembl_id, tf_symbol = source_symbol, layer = level, cytokine = source_category)
  tf_deg <- tf_deg %>% filter(str_detect(layer, "tf-target")) %>% select(-layer) %>% unique()
  
  #join degs together
  if (idx == 1){
    all_tf_deg <- tf_deg
  } else {
    all_tf_deg <- rbind(all_tf_deg, tf_deg)
  }
}

rm(network, network2, tfs, degs, tf_deg)

#### REFORMAT TF DATA ####

# Collapse network by cytokine
all_tfs2 <- all_tfs %>% 
  group_by(tf_ensembl_id, tf_symbol) %>% 
  mutate(target_category = paste0(cytokine, collapse=",")) %>%
  distinct()

# Remove node columns and collapse again
all_tfs2 <- all_tfs2 %>%
  ungroup() %>%
  select (-c(tf_ensembl_id, tf_symbol)) %>%
  group_by(cytokine, target_category) %>%
  summarize(start_count=n())

# Count the number of each target category and join to df
target_count <- all_tfs2 %>%
  ungroup() %>%
  group_by(target_category) %>%
  summarise(target_cat_count = n())
all_tfs3 <- left_join(all_tfs2, target_count)

# Add number column for the target width
all_tfs3 <- all_tfs3 %>% mutate(end_count = start_count/target_cat_count) %>% select(-c(target_cat_count))

# Count the number of each start width and join to df
start_count <- all_tfs3 %>% select(cytokine, start_count) %>%
  ungroup() %>%
  group_by(cytokine) %>%
  summarise(cyto_cat_count = sum(start_count))
all_tfs3 <- left_join(all_tfs3, start_count)

# Edit the start width
all_tfs3 <- all_tfs3 %>% mutate(start_count_n = (10/cyto_cat_count)*start_count) %>% select(-c(cyto_cat_count, start_count))

# Make the output cytokine capital
all_tfs3 <- all_tfs3 %>% ungroup() %>% mutate(cytokine = toupper(cytokine))

#replace the end category col with letter so the labels don't overlap
all_tfs4 <- all_tfs3 %>% select(cytokine,target_category,start_count_n,end_count)

rm(all_tfs2, all_tfs, target_count)

#### PLOT TF DATA ####

# Get list of categories
list_cyto <- all_tfs4[['cytokine']] %>% unique()
list_cap_cyto <- rev(all_tfs4[['target_category']]) %>% unique()
len1 <- length(list_cyto) -1
len2 <- length(list_cap_cyto)-1

# Colours
colours2 <- colours %>% select(category,colour) %>% na.omit()
colours2 <- colours2 %>% filter(category %in% list_cyto |category %in% list_cap_cyto ) 
colours2_v <- setNames(as.character(colours2$colour), colours2$category)

#circos.clear()
#circos.par(start.degree =166, gap.after=c(rep(0, 4), 20, rep(0, 10), 20))
#chordDiagram(all_tfs4, big.gap=20, grid.col=colours2_v, annotationTrack = c("name", "grid"), 
#             annotationTrackHeight = c(0.03, 0.09), order = c("IFNG","IL13","IL17","IL22","TNFA","k", "j","i","h","g","f","e","d","c","b","a"))

circos.clear()
circos.par(start.degree =153, gap.after=c(rep(0, len1), 20, rep(0, len2), 20))
chordDiagram(all_tfs4, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.07),
             big.gap=20, grid.col=colours2_v,annotationTrackHeight = c(0.07, 1), order = c(list_cyto, list_cap_cyto))

#for(si in get.all.sector.index()) {
#  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
#              facing = "inside", niceFacing = TRUE, col = "black",cex=1.5)
#}


#### REFORMAT DEG DATA ####

# Collapse network by cytokine
all_degs2 <- all_degs %>% 
  group_by(deg_ensembl_id, deg_symbol) %>% 
  mutate(target_category = paste0(cytokine, collapse=",")) %>%
  distinct()

# Remove node columns and collapse again
all_degs2 <- all_degs2 %>%
  ungroup() %>%
  select (-c(deg_ensembl_id, deg_symbol)) %>%
  group_by(cytokine, target_category) %>%
  summarize(start_count=n())

# Count the number of each target category and join to df
target_count <- all_degs2 %>%
  ungroup() %>%
  group_by(target_category) %>%
  summarise(target_cat_count = n())
all_degs3 <- left_join(all_degs2, target_count)

# Add number column for the target width
all_degs3 <- all_degs3 %>% mutate(end_count = start_count/target_cat_count) %>% select(-c(target_cat_count))

# Count the number of each start width and join to df
start_count <- all_degs3 %>% select(cytokine, start_count) %>%
  ungroup() %>%
  group_by(cytokine) %>%
  summarise(cyto_cat_count = sum(start_count))
all_degs3 <- left_join(all_degs3, start_count)

# Edit the start width
all_degs3 <- all_degs3 %>% mutate(start_count_n = (75/cyto_cat_count)*start_count) %>% select(-c(cyto_cat_count, start_count))

# Make the output cytokine capital
all_degs3 <- all_degs3 %>% ungroup() %>% mutate(cytokine = toupper(cytokine))

#replace the end category col with letter so the labels don't overlap
all_degs4 <- all_degs3 %>% select(c(cytokine,target_category,start_count_n,end_count))

rm(all_degs2, target_count)

#### PLOT DEG DATA ####

# Get list of categories
list_cyto <- all_degs4[['cytokine']] %>% unique()
list_cap_cyto <- rev(all_degs4[['target_category']]) %>% unique()
len1 <- length(list_cyto) -1
len2 <- length(list_cap_cyto)-1

# Colours
colours2 <- colours %>% select(category,colour) %>% na.omit()
colours2 <- colours2 %>% filter(category %in% list_cyto |category %in% list_cap_cyto ) 
colours2_v <- setNames(as.character(colours2$colour), colours2$category)

#circos.clear()
#circos.par(start.degree =166, gap.after=c(rep(0, 4), 20, rep(0, 10), 20))
#chordDiagram(all_tfs4, big.gap=20, grid.col=colours2_v, annotationTrack = c("name", "grid"), 
#             annotationTrackHeight = c(0.03, 0.09), order = c("IFNG","IL13","IL17","IL22","TNFA","k", "j","i","h","g","f","e","d","c","b","a"))

circos.clear()
circos.par(start.degree = 108, gap.after=c(rep(0, len1), 20, rep(0, len2), 20))
chordDiagram(all_degs4, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.07),
             big.gap=20, grid.col=colours2_v,annotationTrackHeight = c(0.07, 1), order = c(list_cyto, list_cap_cyto))
#for(si in get.all.sector.index()) {
#  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
#              facing = "inside", niceFacing = TRUE, col = "black",cex=1.2)
#}



#### REFORMAT TF-DEG DATA ####

# Collapse tf network by cytokine
all_tf_only <- all_tf_deg %>% 
  select(-c(deg_ensembl_id,deg_symbol)) %>%
  distinct() %>%
  group_by(tf_ensembl_id, tf_symbol) %>% 
  mutate(tf_category = paste0(cytokine, collapse=",")) %>%
  distinct()

# Collapse deg network by cytokine
all_deg_only <- all_tf_deg %>% 
  select(-c(tf_ensembl_id,tf_symbol)) %>%
  distinct() %>%
  group_by(deg_ensembl_id,deg_symbol) %>% 
  mutate(deg_category = paste0(cytokine, collapse=",")) %>%
  distinct()

# Join the tf and deg cytokine categories to the table
all_tf_deg2 <- left_join(all_tf_deg, all_deg_only)
all_tf_deg2 <- left_join(all_tf_deg2, all_tf_only)

# Get count of each tf-deg interaction
all_tf_deg3 <- all_tf_deg2 %>% select(tf_category, deg_category) %>%
  group_by(tf_category, deg_category) %>%
  summarize(tf_deg_count=n())


# Get total number of TFs in each category
num_tfs <- all_tf_only %>% group_by(tf_category) %>% summarise(num_tfs_cat = n())
# Get total size of tf category boxes and append to df
total = sum(num_tfs$num_tfs_cat)
num_tfs <- num_tfs %>% mutate(total_tf_cat = (num_tfs_cat/total)*100)
all_tf_deg3 <- left_join(all_tf_deg3, num_tfs)

# Get tf-deg outgoing size proportion
tf_count <- all_tf_deg3 %>%
    ungroup() %>%
    group_by(tf_category) %>%
    summarise(tf_total = sum(tf_deg_count))
all_tf_deg4 <- left_join(all_tf_deg3, tf_count)
all_tf_deg4 <- all_tf_deg4 %>% mutate("tf_deg_count/tf_total" = tf_deg_count/tf_total)

#Get total size of outgoing
all_tf_deg4 <- all_tf_deg4 %>% mutate(total_outgoing_size = `tf_deg_count/tf_total`*total_tf_cat)


# Get total number of DEGs in each category
num_degs <- all_deg_only %>% group_by(deg_category) %>% summarise(num_degs_cat = n())
# Get total size of deg category boxes and append to df
total2 = sum(num_degs$num_degs_cat)
num_degs <- num_degs %>% mutate(total_deg_cat = (num_degs_cat/total2)*100)
all_tf_deg5 <- left_join(all_tf_deg4, num_degs)

# Get tf-deg incoming size proportion
deg_count <- all_tf_deg5 %>%
  ungroup() %>%
  group_by(deg_category) %>%
  summarise(deg_total = sum(tf_deg_count))
all_tf_deg5 <- left_join(all_tf_deg5, deg_count)
all_tf_deg5 <- all_tf_deg5 %>% mutate("tf_deg_count/deg_total" = tf_deg_count/deg_total)

#Get total size of incoming
all_tf_deg5 <- all_tf_deg5 %>% mutate(total_incoming_size = `tf_deg_count/deg_total`*total_deg_cat)

#Select final column
all_tf_deg6 <- all_tf_deg5 %>% select(tf_category, deg_category, total_outgoing_size, total_incoming_size)
#convert deg category to capitals so different bars in the plot
all_tf_deg6$deg_category <- sapply(all_tf_deg6$deg_category, toupper)

#### PLOT TF DATA ####

# Get list of categories
list_cyto <- all_tf_deg6[['tf_category']] %>% unique()
list_cap_cyto <- rev(all_tf_deg6[['deg_category']]) %>% unique()
len1 <- length(list_cyto) -1
len2 <- length(list_cap_cyto)-1

# Colours
colours2 <- colours %>% select(category,colour) %>% na.omit()
colours2 <- colours2 %>% filter(category %in% list_cyto |category %in% list_cap_cyto ) 
colours2_v <- setNames(as.character(colours2$colour), colours2$category)


circos.clear()
circos.par(start.degree = 170, gap.after=c(rep(0, len1), 20, rep(0, len2), 20))
chordDiagram(all_tf_deg6, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.07),
             big.gap=50, grid.col=colours2_v,annotationTrackHeight = c(0.07, 1),
             order = c(list_cyto, list_cap_cyto))
               