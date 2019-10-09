# Script to visualise the cytokine networks as one circos plot with each circle representing one layer of the networks.
#
# Input: 
#     working_dir - Must contain all cytokine-deg networks to be visualised. Network file names must end with 'omnip_deg_network.txt'.
#     symbol_convert - Biomart downloaded conversion file from human Ensembl IDs to gene symbols.
#     colours_file - Text file containing the colour (hex code) for each cytokine group eg. 'Il13,tnfa'. Group in column 'category' and hex code in column 'colour'.
# Output:
#     The script creates one circos plot but it it doesn't save it automatically.

#### SETUP ####

library(circlize)
library(tidyverse)
library(stringr)

setwd("../input_data/causal_networks/cytokine_deg_networks/")
symbol_convert = "../input_data/id_conversion/mart_export_ensembl_symbol_human_100519.txt" 
colours_file = "../input_data/plot_colours_all.txt"

#### REFORMAT NETWORKS ####

# Import data
all_files <- list.files(pattern = "omnip_deg_network.txt$")
colours = read.delim(colours_file)

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
  
  #Get rows as nodes not interactions
  sources <- network2 %>% select(ensembl_id = source_ensembl_id, symbol = source_symbol, layer = source_layer, category = source_category)
  targets <- network2 %>% select(ensembl_id = target_ensembl_id, symbol = target_symbol, layer = target_layer, category =target_category)
  network3 <- rbind(sources, targets)
  network3 <- distinct(network3)
  
  #join dfs together
  if (idx == 1){
    all_nets <- network3
  } else {
    all_nets <- rbind(all_nets, network3)
  }
    
}

rm(network, network2, network3, sources, targets)

# Collapse network by cytokine
all_nets2 <- all_nets %>% 
  group_by(ensembl_id, symbol, layer) %>% 
  mutate(category = paste0(category, collapse=",")) %>%
  distinct()

# Add missing gene symbols
symbols <- read.delim(symbol_convert, sep = "\t",header=TRUE,na.strings=c("","NA"))
#symbols <-  symbols %>% group_by(Gene.stable.ID) %>% 
#  mutate(Gene.name = paste0(Gene.name, collapse=",")) %>%
#  distinct()
all_nets3 <- left_join(all_nets2, symbols, by = c("ensembl_id" = "Gene.stable.ID"))
all_nets3 <- all_nets3 %>% ungroup() %>% select(-symbol)

# Regroup by symbol
all_nets4 <- all_nets3 %>% 
  group_by(Gene.name, category, layer) %>%
  mutate(ensembl_id = paste0(ensembl_id, collapse=",")) %>%
  distinct()

# Join the colours column
all_nets5 <- left_join(all_nets4, colours, by = c("category" = "category"))

rm(all_nets2, all_nets, all_nets3, all_nets4, symbols)

#### PLOT ####

# Subtables by type
receptors <- all_nets5 %>% filter(layer == "receptor") %>% ungroup() %>%
  select(-Gene.name, -layer, -ensembl_id) %>%
  arrange(category) %>%
  group_by(category, colour)  %>%
  summarize(count=n())
cytokines <- all_nets5 %>% filter(layer == "cytokine")
signals <- all_nets5 %>% 
  filter(layer == "signalling") %>%
  ungroup() %>%
  select(-Gene.name, -layer, -ensembl_id) %>%
  arrange(category) %>%
  group_by(category, colour)  %>%
  summarize(count=n())
tfs <- all_nets5 %>% 
  filter(layer == "tf") %>%
  ungroup() %>%
  select(-Gene.name, -layer, -ensembl_id) %>%
  arrange(category) %>%
  group_by(category, colour)  %>%
  summarize(count=n())
targets <- all_nets5 %>% 
  filter(layer == "target") %>%
  ungroup() %>%
  select(-Gene.name, -layer, -ensembl_id) %>%
  arrange(category) %>%
  group_by(category, colour)  %>%
  summarize(count=n())

circos.clear()

# Circle of targets
len1 <- length(targets$category %>% unique())
circos.par(gap.degree = 0, cell.padding = c(0.02, 0, 0.02, 0), "track.height" = 0.125 ) #No gap between sections
circos.initialize(factors = as.factor(targets$category), xlim = cbind(rep(0,len1), targets$count))
circos.track(ylim=c(0,1), bg.col = as.character(targets$colour)) #,panel.fun = function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing='outside', niceFacing=TRUE, adj = c(0.5, 5))})
circos.clear()

# circle of tfs
par(new = TRUE)
len2 <- length(tfs$category %>% unique())
circos.par("canvas.xlim" = c(-1.25, 1.25), "canvas.ylim" = c(-1.25, 1.25), gap.degree= 0, cell.padding = c(0.02, 0, 0.02, 0), "track.height" = 0.16)
circos.initialize(factors = tfs$category, xlim = cbind(rep(0,len2), tfs$count))
circos.track(ylim = c(0, 1), bg.col = as.character(tfs$colour))
circos.clear()

# circle of signalling
len3 <- length(signals$category %>% unique())
par(new = TRUE)
circos.par("canvas.xlim" = c(-1.7, 1.7), "canvas.ylim" = c(-1.7, 1.7), gap.degree= 0,cell.padding = c(0.02, 0, 0.02, 0), "track.height" = 0.22)
circos.initialize(factors = signals$category, xlim = cbind(rep(0,len3), signals$count))
circos.track(ylim = c(0, 1), bg.col = as.character(signals$colour))
circos.clear()

# circle of receptors
len4 <- length(receptors$category %>% unique())
par(new = TRUE)
circos.par("canvas.xlim" = c(-2.7, 2.7), "canvas.ylim" = c(-2.7, 2.7), gap.degree= 0, cell.padding = c(0.02, 0, 0.02, 0), "track.height" = 0.3)
circos.initialize(factors = receptors$category, xlim = cbind(rep(0,len4), receptors$count))
circos.track(ylim = c(0, 1), bg.col = as.character(receptors$colour))
circos.clear()

# circle of cytokines
par(new = TRUE)
circos.par("canvas.xlim" = c(-5.9, 5.9), "canvas.ylim" = c(-5.9, 5.9), gap.degree= 0, cell.padding = c(0.02, 0, 0.02, 0),  "track.height" = 0.6)
circos.initialize(factors = cytokines$Gene.name, xlim = c(0, 1))
circos.track(ylim = c(0, 1), bg.col = as.character(cytokines$colour))#,panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, facing='outside', niceFacing=TRUE, adj = c(0.5, 0.5), cex=1  )})
circos.clear()

