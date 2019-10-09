"""
Script to filter the input TF-DEG networks for interactions where the TF
passes the given corrected p value filter based on the output of the TF filtering
internal tool (which must be run prior to this).

Input:
    1. Path to the folder containing all output files from the TF filtering internal tool.
    The files must be named '{cytokine}_results.txt'.
    2. Path to the folder containing all the input TF-DEG networks to be filtered. The
    networks must be named '{cytokine}_tf_deg_network'.

Output:
    1. tf_filter_{pcorr}/{cytokine}_tfs_filt_{pcorr}.txt - List of TFs passing the filter.
    2. tf_filter_{pcorr}/cytokine}_tf_deg_network_chatfilt_{pcorr}.txt - filtered TF-DEG network

"""

import pandas as pd
import os

pcorr = 0.1 # (< pcorr)
tffilt_path = "input_data/tf_filtering/tf_filter_results/"
unfilt_net_path = "input_data/causal_networks/tf_deg_networks/no_tf_filtering/"
output_path = f"input_data/causal_networks/tf_deg_networks/"
df_names = ["ifng", "il22", "tnfa", "il13"]


# For each cytokine
for cytokine in df_names:

    # Open unfiltered network
    filename = f"{cytokine}_tf_deg_network.txt"
    filepath = f"{os.path.join(unfilt_net_path, filename)}"
    unfilt_net = pd.read_csv(filepath, sep = "\t")

    # Get list of regualtors in the network
    regs = pd.Series(unfilt_net['source_ensembl_id']).unique()

    # Open chat output
    filename1 = f"{cytokine}_results.txt"
    filepath1 = f"{os.path.join(tffilt_path, filename1)}"
    chat_output = pd.read_csv(filepath1, sep = "\t")

    # Get the chat results for the regulators in the network and filter for p corrected value
    chat_regs = chat_output[(chat_output['Node'].isin(regs)) & (chat_output['p-corr'] < pcorr)]

    # Save these filtered TF lists
    folder = f"tf_filter_{pcorr}"
    dirName = f"{os.path.join(output_path, folder)},"
    if not os.path.exists(dirName):
        os.mkdir(dirName)
    filename2 = f"{cytokine}_tfs_filt_{pcorr}.txt"
    filepath2 = f"{os.path.join(dirName, filename2)}"
    chat_regs.to_csv(filepath2, sep="\t", index=False)

    # Filter the network for only these regulators
    filt_net = unfilt_net[unfilt_net['source_ensembl_id'].isin(chat_regs['Node'])]
    #print(filt_net)

    # Save the fitlered network
    dirName2 = f"{os.path.join(dirName, 'Filtered_networks')},"
    if not os.path.exists(dirName2):
        os.mkdir(dirName2)
    filename3 = f"{cytokine}_tf_deg_network_filt_{pcorr}.txt"
    filepath3 = f"{os.path.join(dirName2, filename3)}"
    filt_net.to_csv(filepath3, sep="\t", index=False)
    
