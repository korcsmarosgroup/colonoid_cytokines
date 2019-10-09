"""
Extract paths to shared tfs and combine into one network with nodes labelled according to which networks they are in.

Input:
    network_dir: Path to cytokine-tf networks.
    sharedtfs_path: Path to text file with the cytokine-tf direct connections where the tfs are present in >1 network.
            The cytokine-tf direct file is generated in the script 'cytokine_to_tf_or_deg_direct.py' and then the
            cytokine-tf interactions only occuring in one network (one cytokine) are removed manually.
Output:
    Network file of all interactions between all cytokines and the shared TFs. Each interaction is labelled with which
    network it occurs in.
"""
###############
# Setup
###############

# import
import pandas as pd
import os
import networkx as nx

# input files
network_dir = "../input_data/causal_networks/cytokine_tf_networks/"
sharedtfs_path = "../input_data/causal_networks/cytokine_tf_networks/shared_tfs/AllCytokines_to_tf_shared_network.txt"
out_path = "../input_data/causal_networks/cytokine_tf_networks/shared_tfs/Network_all_paths_to_shared_tfs.txt"

# Load tfs
tfs = pd.read_csv(sharedtfs_path, sep="\t")
#tf_list = tfs['tf_ensembl_id'].tolist()

########################
# Get network filepaths
########################

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

########################
# Get signalling paths
########################

# Empty dataframe for the edges
paths_net_all = pd.DataFrame()

# Iterate filepaths to get one large network
for idx, filename in enumerate(netfolder.files):
    print(filename)
    # Open the networks
    filepath = f"{os.path.join(network_dir, filename)}"
    network_df = pd.read_csv(filepath, sep="\t")

    # Convert to network
    net = nx.from_pandas_edgelist(network_df, 'source_ensembl_id', 'target_ensembl_id', create_using=nx.DiGraph())

    # Get cytokine of current network
    cytokine_dict = {
        "ifng": "ENSG00000111537",
        "il17": "ENSG00000112115",
        "tnfa": "ENSG00000228321,ENSG00000206439,ENSG00000204490,ENSG00000232810,ENSG00000230108,ENSG00000223952",
        "il13": "ENSG00000169194"}
    cyto = filename.split("_")[0]

    # Get TFs targeted by specific cytokine
    cyto_tfs_df = tfs[tfs['cytokine_all'].str.contains(cyto)]
    tf_list = cyto_tfs_df['tf_ensembl_id'].tolist()

    # Get paths to tfs
    if (cyto == "tnfa"):
        # Account for having >1 id
        value1 = cytokine_dict[cyto].split(",")
        for val in value1:
            print("\n" + val + "\n")
            for tf in tf_list:
                print("\n" + tf + "\n")
                # Get all the shortest paths between source and target
                try:
                    for p in nx.all_simple_paths(net, source=val, target=tf):

                            # Convert path to edges and append to df
                            for idx, item in enumerate(p):

                                if (len(p) - 1) == idx:
                                    continue
                                next = idx + 1
                                edges = {'source': [item],
                                         'target': [p[next]],
                                         'cyto': [cyto]}
                                df_edges = pd.DataFrame(edges, columns=['source', 'target', 'cyto'])
                                # Append the edges to the others for the cytokine
                                paths_net_all = paths_net_all.append(df_edges)
                except nx.NetworkXNoPath:
                    print("no paths found for " + cyto + ", source= " + val + ", target = " + tf)
                    continue
    else:
        for tf in tf_list:
            print("\n"+tf+"\n")
            # Get all the shortest paths between source and target
            try:
                for p in nx.all_simple_paths(net, source=cytokine_dict[cyto], target=tf):

                    # Convert path to edges and append to df
                    for idx, item in enumerate(p):

                        if (len(p) - 1) == idx:
                            continue
                        next = idx + 1
                        edges = {'source': [item],
                                 'target': [p[next]],
                                 'cyto': [cyto]}
                        df_edges = pd.DataFrame(edges, columns=['source', 'target', 'cyto'])
                        # Append the edges to the others for the cytokine
                        paths_net_all = paths_net_all.append(df_edges)
            except nx.NetworkXNoPath:
                print("no paths found for " + cyto + ", source= " + val + ", target = " + tf)
                continue

##################################
# Clean and save signalling paths
##################################

# Remove duplicates in paths_net_all
if not paths_net_all.empty:
    paths_net_all = paths_net_all.drop_duplicates()

    # Collapse by cytokine network
    paths_net_all_group = paths_net_all.groupby(['source', 'target'], as_index=False).agg(
            {'cyto': ','.join}).reset_index()

    paths_net_all_group.to_csv(out_path, sep='\t', index=False, header=True)