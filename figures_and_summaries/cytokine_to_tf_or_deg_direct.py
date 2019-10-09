"""
Script to extract the direct links between cytokine and tfs or DEGs from the causal networks and to specify
if the cytokines share tfs.

Input:
    network_dir: Path to folder containing the cytokine-DEG networks output from the 'network_generation.py' script.
                 Network files must be named '{cytokine}_omnip_deg_network.txt'.
    type: Either 'tf' or 'deg' to specify whether to generate cytokine-deg or tf direct links.
Output:
    Tab delimited text file with the TFs or DEGs and which cytokine networks they are in.

"""

import os
import pandas as pd

# Parameters and inputs
network_dir = "../input_data/causal_networks/cytokine_deg_networks/"
output_path =  "AllCytokines_to_tf_direct_network.txt"
type = "tf" #"tf" or "deg"

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
tf_dict = {}

# Iterate filepaths
for idx, filename in enumerate(netfolder.files):

    # Open the network
    filepath = f"{os.path.join(network_dir, filename)}"
    network = pd.read_csv(filepath, sep="\t")

    # Get cytokine name
    cytokine = filename.replace("_omnip_deg_network.txt", "")

    # Get all the tf-target edges
    tftgs = network[network['level']=='tf-target']

    if type == "tf":
        tftgs_dict = tftgs.set_index('source_ensembl_id')['source_symbol'].to_dict()
    else:
        tftgs_dict = tftgs.set_index('target_ensembl_id')['target_symbol'].to_dict()

    # Add cytokine to value in dict
    if idx == 0:
        for key, value in tftgs_dict.items():
            tf_dict[key] = f"{value}:{cytokine}"
    else:
        for key, value in tftgs_dict.items():
            if key in tf_dict:
                curr_val = tf_dict[key]
                new_val = f"{curr_val},{cytokine}"
                tf_dict[key] = new_val
            else:
                tf_dict[key] = f"{value}:{cytokine}"
                
# Convert the dict as a dataframe
tf_df = pd.DataFrame.from_dict(tf_dict, orient='index')

# Reformat the df
tf_df = tf_df.reset_index()
tf_df2 = pd.DataFrame(tf_df[0].str.split(':',1).tolist(),columns = ['symbol','cytokine'])
tf_df3 = pd.concat([tf_df, tf_df2], axis=1).drop([0], axis=1)
if type == "tf":
    tf_df3.columns = ['tf_ensembl_id', 'tf_symbol','cytokine_all']
else:
    tf_df3.columns = ['deg_ensembl_id', 'deg_symbol', 'cytokine_all']

# One line per cytokine so can add to cytoskape as network
tf_df4 = tf_df3['cytokine_all'].str.split(',', expand=True).stack().reset_index(level=0).set_index(
    'level_0').rename(columns={0: 'cytokine'}).join(tf_df3) #.drop('index', 1))
tf_df4 = tf_df4.reset_index(level=0, drop=True)

# Add column for cytokine symbol - to enable cytoskape visualisation
tf_df4['cytokine_symbol'] = tf_df4['cytokine'].str.upper()

# Save output
tf_df4.to_csv(output_path, sep = "\t", index=False)