# Count number of DEGs in each network
#
# Input:
#     folder: Path to folder containing the cytokine-DEG networks output from the 'network_generation.py' script.
#           Network files must be named '{cytokine}_omnip_deg_network.txt'
# Output:
#     Script print to screen the number of DEGs in each cytokine networks.

library(dplyr)

# Input
folder <- "../input_data/causal_networks/cytokine_deg_networks/"
file <- "_omnip_deg_network.txt"
cytokines <- c("ifng","il13","il17","tnfa")

# Iterate cytokines
for (cyto in cytokines){
  filepath <- paste(folder, cyto, file, sep="")
  cyto_df <- read.csv(filepath, sep="\t", header= TRUE)
  # Get distinct targets of tf-target interactions
  cyto_df2 <- cyto_df %>% filter(level == "tf-target") %>% select(target_ensembl_id) %>% distinct()
  # Print the cytokine and the number of DEGs
  print (paste(cyto, ": ", nrow(cyto_df2)))
}
