# Script to load all networks of shared interactions (one for each combinationof cytokines) into Cytoscape,
# automatically do the visualisation and save as pdf images
# NB. You must have Cytoscape open for this script to work.

# Input:
#       network - All interactions which are shared between cytokine networks with the 'Cytokines' column indicating which networks they are present in. 
#                 Created in the combine_network_data.py script (and edited to remove interactions in only one network).
#       colours - Text file containing the colour (hex code) for each cytokine group eg. 'Il13,tnfa'. Group in column 'category' and hex code in column 'colour'.
#       symbols - Biomart downloaded conversion file from human Ensembl IDs to gene symbols.
# Output:
#       Saves a PDF of each network to file (each network has all interactions shared between one combination of cytokines eg. tnfa and ifng, il13 and ifng.

##### Set up #####

library(dplyr)
library(RCy3)
cytoscapePing()

# Input paths
network <- read.csv("../input_data/causal_networks/cytokine_deg_networks/shared_interactions/Interactions_2ormore.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
colours <- read.csv("../input_data/plot_colours_all.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
symbols <- read.csv("../input_data/id_conversion/mart_export_ensembl_symbol_human_100519.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

# Get node network membership
network_s <- network %>% dplyr::select(c(node = source, Cytokines)) %>% distinct()
network_t <- network %>% dplyr::select(c(node = target, Cytokines)) %>% distinct()
network_nodes <- rbind(network_s, network_t) %>% distinct()
network_nodes <- network_nodes %>% group_by(node) %>% 
  mutate(node_membership = paste0(Cytokines, collapse = "|")) %>% 
  dplyr::select(c(node, node_membership)) %>% distinct()
network_nodes <- as.data.frame(network_nodes)

##### Read in network #####

# Create network in Cytoscape
# Ensure there are no files already in cytoscape
cyto_net<-createNetworkFromDataFrames(edges=network, title="shared_ints_all", collection="Cytokine Networks")

# Add gene symbols
loadTableData(symbols, data.key.column = "Gene.stable.ID", table.key.column = "shared name", network = "shared_ints_all")

# Add node membership
loadTableData(network_nodes, data.key.column = "node", table.key.column = "shared name", network = "shared_ints_all")

##### Create subnetworks #####

# Categories
cats <- network$Cytokines %>% unique()

# Create subnetworks
for (item in cats){
  createColumnFilter("cytokine_filter", 'Cytokines', item, "IS", network ="shared_ints_all", type = "edges", anyMatch =FALSE)
  #selectNodesConnectedBySelectedEdges()
  createSubnetwork(subnetwork.name = item, edges = "selected")
}

##### Visual parameters #####

# Node defaults
setNodeColorDefault('#FFFFFF') # set to grey
setNodeShapeDefault('ELLIPSE')
setNodeSizeDefault(60)
lockNodeDimensions(TRUE)
setNodeLabelDefault(new.label = "")
setNodeBorderColorDefault('#000000')
setNodeBorderWidthDefault("10")
setNodeFontFaceDefault(new.font="AlBayan-Bold")

for (item in cats){
  col = colours %>% filter(category == item)
  # Layouts
  layoutNetwork(layout.name ="kamada-kawai", network = item)
  # Node colour
  list <- getAllNodes(network = item)
  setNodePropertyBypass(list, col$colour, "NODE_FILL_COLOR", network=item)
}

# MANUALLY CHECK AND FINE TUNE

##### Export images #####

for (item in cats){
  fitContent(network=item)
  exportImage(filename = paste0(item, "_shared_ints.pdf", sep = ""), type="PDF", network=item)
}
