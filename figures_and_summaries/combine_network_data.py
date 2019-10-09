"""
Script to create table of which nodes/interactions are shared between which networks. It outputs a table of
nodes/interactions as well as a venn diagram of the overlaps.

Input:
    network_dir: Path to cytokine-tf networks.
    type_id:  'nodes' or 'interactions' depending on whether you want to summarise nodes or interactions of the networks

Output:
    - Table of each node/interaction and columns indicting whether that node/interaction is present in each network -
     and a count of how many networks it is in.
    - A Venn diagram of how many nodes/interactions are shared between each of the networks.

"""

#########
# Setup
#########

import os
import pandas as pd
from venn import venn

network_dir = "../input_data/causal_networks/cytokine_deg_networks/"
output_path = "../input_data/causal_networks/cytokine_deg_networks/Network_summary_nodes.txt"
venn_out_path = "../input_data/causal_networks/cytokine_deg_networks/Venn_summary_nodes.svg"
type_id = "nodes" # "interactions" or "nodes"

####################
# Get network files
####################

class FolderHandler:
    def __init__(self, path):

        if path[:-1] == "/":
            self.path = path[0:-1]
        else:
            self.path = path

        self.files1 = os.listdir(self.path)

        self.files = []
        for names in self.files1:
            if names.endswith("network.txt"):
                self.files.append(names)

# Instantiate class netfolder. init functions called automatically
netfolder = FolderHandler(network_dir)
list_dict = {}

###################
# Process networks
###################

# Iterate filepaths
for idx, filename in enumerate(netfolder.files):

    # Open the networks
    filepath = f"{os.path.join(network_dir, filename)}"
    network = pd.read_csv(filepath, sep="\t")

    if type_id == "nodes":
        # get source nodes
        output_df = pd.DataFrame(network['source_ensembl_id']).drop_duplicates().rename(
            columns={'source_ensembl_id': type_id})
        # get target nodes
        targets = pd.DataFrame(network['target_ensembl_id']).drop_duplicates().rename(
            columns={'target_ensembl_id': type_id})
        # addend source and targets
        output_df = output_df.append(targets).drop_duplicates()
    else:
        # Combine the source and targets nodes as a path
        output_df = pd.DataFrame(network[['source_ensembl_id','target_ensembl_id']])
        # Concatenate source and target
        output_df[type_id] = output_df["source_ensembl_id"].map(str) + ":" + output_df["target_ensembl_id"].map(str)
        output_df  = output_df.drop(['source_ensembl_id','target_ensembl_id'], axis=1)

    # Convert to set via list
    output_df_list = output_df[type_id].tolist()
    output_df_set = set(output_df_list)
    # Add set to dictionary for venn. key=name, value =list
    cytokine = filename.replace("_omnip_deg_network.txt","")
    #print(cytokine + " : " + len(output_df_set))
    l = len(output_df_set)
    print(l)
    list_dict[cytokine] = output_df_set
    # add columns of 1's for first network
    output_df[filename] = 1

    if idx == 0:
        # For the first network
        output_df_main = output_df
    else:
        # For the other networks merge with main network
        output_df_main = output_df_main.merge(output_df, how='outer')

####################
# Save matrix file
####################

#check no duplicates, convert nan to 0
output_df_main = output_df_main.drop_duplicates()
# Put nodes as index column
output_df_main = output_df_main.set_index(type_id)
#  Add sum columns
output_df_main['count'] = output_df_main.count(axis='columns')
# Sort by count
output_df_main = output_df_main.sort_values('count', ascending=False)

# Save as file
output_df_main.to_csv(output_path, sep = "\t", na_rep='Na')

##############
# Create Venn
##############

# Colours specification
cmap = ["#FFFF00","#5C85FF", "#FF0033","#6d6c6c"]

plot = venn(list_dict, cmap=cmap)
fig = plot.get_figure()
fig.savefig(venn_out_path,format='svg', dpi=1200)