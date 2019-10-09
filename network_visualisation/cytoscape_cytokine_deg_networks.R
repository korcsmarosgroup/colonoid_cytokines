# Script to load networks into Cytoscape and automatically do the visualisation
# NB. You must have Cytoscape open for this script to work.
#
# Input:
#     network: Path to folder containing the cytokine-DEG networks output from the 'network_generation.py' script.
#           Network files must be named '{cytokine}_omnip_deg_network.txt'
#     degs: Path to folder containing the processed differentially expressed gene lists for each cytokine dataset.
#           These files are created and saved during the 'network_generation.py' script. 
#           Files must be named '{cytokine}_unfiltered.txt'.
#
# Output: 
#     The script loads the networks into Cytoscape but does not save them as Cys files or as images, this must be done manually.

##### Set up #####

library(dplyr)
library(RCy3)
cytoscapePing()

# Input paths
network <- "../input_data/causal_networks/cytokine_deg_networks/"
degs <- "../differential_expression_datasets/colonoids_processed/"

# Networks to process
cytokines <- c("ifng", "il17", "il13","tnfa")

##### Read in networks #####

for (cyto in cytokines){
  
  # Filepaths
  net_path <- paste(network, cyto, "_omnip_deg_network.txt", sep ="")
  degs_path <-  paste(degs, cyto, "_unfiltered.txt", sep ="")
  
  # Import network as data frame
  cyto_df <- read.csv(net_path, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get table of edges
  cyto_edges <- data.frame(source=cyto_df$source_ensembl_id,
                      target=cyto_df$target_ensembl_id,
                      level = sapply(strsplit(cyto_df$level,","), `[`, 1),
                      stringsAsFactors=FALSE)
  
  # Create network in Cytoscape
  cyto_net<-createNetworkFromDataFrames(edges=cyto_edges, title=paste(cyto, " Network", sep = ""), collection="Cytokine Networks")
  
  
  # Get the layer position of each node (highest possible position in the signalling pathway)
  cyto_edges_f <- cyto_edges %>% select(c(-source)) %>% group_by(node=target) %>% summarise(levels=paste(level, collapse=","))
  cyto_edges_f <- cyto_edges_f %>% mutate(node_layer = ifelse(grepl('cytokine-receptor',levels), "receptor", 
                                                                  ifelse(grepl('receptor-interm1', levels), "signalling",
                                                                         ifelse(grepl('interm-interm', levels), "signalling",
                                                                                ifelse(grepl('final_interm-tf', levels), "tf",
                                                                                       "deg")))))
  to_add <- subset(cyto_edges, !(source %in% cyto_edges_f$node)) %>% select(c(node =source)) %>% mutate(levels = "cytokine", node_layer = "cytokine")
  cyto_edges_f <- rbind(cyto_edges_f, to_add) %>% ungroup()
  cyto_edges_f2 <- as.data.frame(cyto_edges_f)
  loadTableData(cyto_edges_f2, data.key.column = "node", table.key.column = "shared name", network = paste(cyto, " Network", sep = ""), table = "node")
  
  # Add lfc, q val and symbols
  cyto_degs_df <- read.csv(degs_path, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  loadTableData(cyto_degs_df, data.key.column = "ENSEMBL.ID", table.key.column = "shared name", network = paste(cyto, " Network", sep = ""))
  
}

##### Visual parameters #####

# Node defaults
deleteStyleMapping("NODE_LABEL")
setNodeColorDefault('#FFFFFF') # set to grey
setNodeShapeDefault('ELLIPSE')
setNodeSizeDefault(20)
lockNodeDimensions(TRUE)
setNodeLabelDefault(new.label = "")
setNodeBorderColorDefault('#000000')
setNodeBorderWidthDefault("6")

# Edge defaults
setEdgeTargetArrowShapeDefault("Delta")
setEdgeLineWidthDefault("4")
matchArrowColorToEdge(TRUE)

# Node mappings
setNodeColorMapping("adj.P.Val", c(0.01,0.010000001), c('#000000','#FFFFFF'), mapping.type = "c")
setNodeSizeMapping("logFC", c(0, 7), c(20,100), mapping.type = "c")

# Edge mapping
setEdgeColorMapping("level", c("cytokine-receptor","receptor-interm1","interm-interm","final_interm-tf","tf-target"),
                    c("#33a02c","#b2df8a","#1f78b4", "#1f78b4","#a6cee3"), mapping.type = "d")

# Layouts
setLayoutProperties('force-directed',list(defaultSpringLength=50, defaultSpringCoefficient=6E-04))

