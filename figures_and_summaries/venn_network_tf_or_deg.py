"""
Create venn diagram of the tfs or degs targeted by each cytokine.

Input:
    network: Path to cytokine-tf or cytokine-deg direct network generated in the script Cytokine_to_tf_or_deg_direct.py.
    type_data: 'tf' or 'deg' to plot the tf or deg Venns.

Output:
    Venn diagram of the overlaps of DEGs or TFs targeted by each cytokine.
"""
import os
import pandas as pd
from venn import venn

# Inputs and parameters
network = "../input_data/causal_networks/cytokine_deg_networks/cytokine_deg_direct/AllCytokines_to_tf_direct_network.txt"
venn_out_path = "Venn_cytokines_to_tf_nodes.svg"
type_data = "tf" # 'tf' or 'deg'

# Load the network
network = pd.read_csv(network, sep="\t")

if type_data == 'tf':
    # Get cytokine and target column only
    network2 = network.drop(['tf_symbol', 'cytokine_all', 'cytokine_symbol'], axis=1).drop_duplicates()
    # Group by cytokine and list the targets as a set
    network2 = network2.groupby('cytokine')['tf_ensembl_id'].apply(set)
else:
    # Get cytokine and target column only
    network2 = network.drop(['deg_symbol', 'cytokine_all', 'cytokine_symbol'], axis = 1).drop_duplicates()
    # Group by cytokine and list the targets as a set
    network2 = network2.groupby('cytokine')['deg_ensembl_id'].apply(set)

# Convert to dictionary
net_dict = network2.to_dict()
#print(network2['ifng'])

# Create Venn
# Colours specification
cmap = ["#FF0033","#FFFF00","#6d6c6c","#5C85FF"]
# Save plot
plot = venn(net_dict, cmap=cmap)
fig = plot.get_figure()
fig.savefig(venn_out_path,format='svg', dpi=1200)