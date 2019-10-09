""" Script to connect diff exp genes (from cytokine treated organoids) with the cytokines. Using TFTG network
and omnipath. Run from scripts directory.

If transcription factor filtering is top be used, the script must be run twice with different input settings.
First the script is run to create TF-DEG networks
    - do_all_load, do_tfnet, do_omnipath_filter_exp, do_filter_hpa = True
    - do_tf_filt, do_tf_filt_chat, do_omnipath = False
They are then filtered externally and read into the script upon the second run through,
to add the upstream signalling pathways.
    - do_all_load, do_tf_filt_chat, do_omnipath, do_omnipath_filter_exp, do_filter_hpa = True
    - do_tf_filt, do_tfnet = False

NB. 'chat' refers to the internal tool which filters TF for their influence in the network using the transcriptomics
 data p values and the degree of the nodes. It is based on the cytoscape app 'CHAT'.

Inputs:
    1. Txt files from BioMart of human Uniprot - Ensembl and Ensembl - Gene symbol ID conversions
    2. DoRothEA TF-TG data in TSV format (does not need to be converted to Ensembl IDs first)
    3. OmniPath data in TSV format (prefiltered for correct cytokine receptor interactions)
    4. Excel files containing the differentially expressed genes for the cytokine treated colonoids (from Powell lab)
    5. Excel files containing the expression data for the genes in cytokine treated colonoids (from Powell lab)
    6. Directory path containing TF-DEG networks prefiltered to include only regulators passing the p corrected value
    required for the CHAT based regulator filtering (and for human protein atlas if applied).
    7. TSV file of genes expressed in the colon downloaded form the Human protein atlas.

Outputs:
    1/2. Tab delimited text files of the biomart ID conversion files with 'na's removed.
    3. Tab delimited text files of all the DEGs passing the lfc and adj p cut offs (in the function 'load_degs'.
    4. Tab delimited text file containing the lists of expressed genes based on given fpkm cut off
        (in the 'load_expressed' function).
    5. TF-TG network with Ensembl IDs added
    6. TF-DEG network filtered by input data saved as text file
    7. Cytokine-TF network filtered by input data saved as text file
    8. Cytokine-DEG network filtered by input data saved as text file

"""

##########
# Imports
##########
import networkx as nx
import pandas as pd
import os

########################
# Parameters / cut offs
########################

adjp = 0.01 # adj p value filter for degs (using <)
chat_pcorr = 0.1 # The p corr cut off that has been applied for the tf filter (if applied)
lfc = 0 # lfc value filter for degs (0= no filter) (using > modulus)
exp_filt = 0 # Value for expression of gene (using >). At least 2 of 3 replicates need this fpkm value.
tf_filt_num = 50 # Number of targets a TF must have in the TF-deg network - for filtering of TFs based on degree (>=)

# If changing df names make sure lines 196 (end of load_degs function), 245 (load_expressed function),
# and 661 (within main) are changed to match
df_names = ["ifng", "il17", "tnfa", "il13"]

######################################
# Which sections of the script to run
######################################

do_all_load = True #Set to True if you want to go from the raw files, select False if you have already run the load
# functions and want to simply load their preprocessed outputs (lots quicker) (doesnt load omnipath)
# If you change the adjp or lfc filters this needs to be true to update them
do_tf_filt = False # Set to True if want to filter the TFs based on number of targets in network instead of loading
#  TF filtered files
do_tfnet = False # Set to true to generate the tf-degs networks, else false to read them in (the unfiltered versions)
do_tf_filt_chat = True # Set to true if you want to open tf pre-filtered tftg networks instead of doing the tf filt
#  based on degree
do_omnipath = True # Set to true to carry out the loading and filtering of the omnipath data
do_omnipath_filter_exp = True # Set as true to filter all omnipath proteins/gene as expressed in the dataset
do_filter_hpa = True # Set as true to filter the whole networks for genes expressed in the colon according
# to the human protein atlas. This happens separately for the tf-deg and the cytokine-tf networks.

# Checks we don't have interfering settings
if do_tf_filt and  do_tf_filt_chat:
    sys.exit('Error: You must filter the TFs either by number of targets or using an external TF based filter,'
             ' it cant be both.')
if not do_tf_filt and not do_tf_filt_chat and do_omnipath:
    sys.exit('Error: If you are to run Omnipath analysis, you must filter the TFs either by number of targets or '
             'using an external TF filter.')

###################
# Input file paths
###################

idconversion_path = "../input_data/id_conversion/mart_export_uniprot_ensembl_human_110419.txt" #uniprot-ensembl biomart
idconversion_path_sym = "../input_data/id_conversion/mart_export_ensembl_symbol_human_100519.txt" #ensembl-gene symbol biomart
tftg_path = "../input_data/prior_knowledge_networks/tfregulons_database_v01_20180216__ABCD.tsv" #Gene symbols
omnipath_path = "../input_data/prior_knowledge_networks/omnipath_0.7.111_formatted_cytokinereceptors.tsv"
degs1_path = "../input_data/differential_expression_datasets/colonoids_raw/il17_tnfa_ifng_degs.xlsx"
degs2_path = "../input_data/differential_expression_datasets/colonoids_raw/il13_degs.xlsx"
expressed1_path = "../input_data/expression_datasets/il17_tnfa_ifng_fpkm.xlsx"
expressed2_path = "../input_data/expression_datasets/il13_fpkm.xlsx"
chat_path = "../input_data/causal_networks/tf_deg_networks/with_tf_filtering/"
hpa_path = "../input_data/HPA_tissue_detectable_rna_colon_yes.tsv"

####################
# Output file paths
####################

idconversion_out = "../input_data/id_conversion/mart_export_uniprot_ensembl_human_filtered.txt"
idconversion_sym_out = "../input_data/id_conversion/mart_export_ensembl_symbol_human_100519_filtered.txt"
tftgs_out = "../input_data/prior_knowledge_networks/tfregulons_database_v01_20180216__ABCD_ensembl2.txt"
# Output file paths for the deg lists
ifng_degs_out = f"../input_data/differential_expression_datasets/colonoids_processed/ifng_degs_adjp_{adjp}.txt"
il17_degs_out = f"../input_data/differential_expression_datasets/colonoids_processed/il17_degs_adjp_{adjp}.txt"
tnfa_degs_out = f"../input_data/differential_expression_datasets/colonoids_processed/tnfa_degs_adjp_{adjp}.txt"
il13_degs_out = f"../input_data/differential_expression_datasets/colonoids_processed/il13_degs_adjp_{adjp}.txt"
# Output folder for the expressed gene lists
exp_out_dir = f"../input_data/expression_datasets/fpkm_morethan_{exp_filt}_in2reps/"
# Output file paths for the tf-deg networks
deg_tf_output_dir = f"../input_data/causal_networks/tf_deg_networks/deg_adjp_{adjp}_hpa_dorothea/"
deg_tf_filt_output_dir = f"../input_data/causal_networks/tf_deg_networks/{tf_filt_num}_or_more_targets/"
# Output file paths for omnipath networks
if do_tf_filt_chat:
    omnip_tf_out_dir = f"../input_data/causal_networks/cytokine_tf_networks/"
    # Output paths for whole cytokine to deg networks
    omnip_deg_out_dir = f"../input_data/causal_networks/cytokine_deg_networks/"
elif do_tf_filt:
    omnip_tf_out_dir = f"../input_data/causal_networks/cytokine_tf_networks/"
    # Output paths for whole cytokine to deg networks
    omnip_deg_out_dir = f"../input_data/causal_networks/cytokine_deg_networks/"

#############
# Functions
#############

def load_id_conversion(idconversion_path, idconversion_path_sym):
    """
    Load and preprocess Uniprot-Ensembl conversion file.

    :param idconversion_path: Path to the uniprot-ensembl ID conversion table downloaded from Biomart. Columns
                    'Gene stable ID' = Ensembl ID and 'UniProtKB/Swiss-Prot ID' = uniprot ID.
    :param idconversion_path_sym: Path to the gene symbol-ensembl ID conversion table downloaded from Biomart. Columns
                    'Gene stable ID' = Ensembl ID and 'Gene name' = gene symbol.

    :return: Dataframe of uniprot to ensembl conversion - single to single, and similar for ensembl to gene symbol.
    """

    # Load table
    idconversion = pd.read_csv(idconversion_path)
    idconversion_sym = pd.read_csv(idconversion_path_sym, sep = "\t")

    # Remove lines with NA
    idconversion_f = idconversion.dropna()
    idconversion_sym_f = idconversion_sym.dropna()

    # Save as file
    idconversion_f.to_csv(idconversion_out, sep='\t', header=True, index=False)
    idconversion_sym_f.to_csv(idconversion_sym_out, sep='\t', header=True, index=False)

    return idconversion_f, idconversion_sym_f


def load_tftg_prior_knowledge(tftg_path):
    """
    Load and preprocess the TF-TG prior knowledge

    :param tftg_path: path to the tftg prior knowledge network. Columns named 'TF' and 'target'.

    :return: dataframe of tftg interactions with same column names.
    """

    # Load network
    # setting dtype stops it erroring because some values are missing in these columns
    tftgs = pd.read_csv(tftg_path, sep ='\t',dtype={"which_TFbindingMotif":"str", "which_coexpression":"str"})

    return tftgs


def load_degs(degs1_path,degs2_path):
    """
    Load and preprocess the differential expression data.

    :param degs1_path: Path to the excel file with DEG data from IL17, IFNg and TNFa (from Nick and Polychronis)

    :param degs2_path: Path to the excel file with DEG data from IL13 (from Nick and Polychronis)

    :return: list of filtered dataframes of degs filtered based on the adjp and logfc cut offs. lfc in column 'logFC'
            and adj p in column 'adj.P.Val'.

    :output: saves tab delimited text files of all the DEGs passing the lfc and adj p cut offs.
    """

    # Load unfiltered deg data
    data_ifng = pd.read_excel(degs1_path, sheet_name='DEAnalysisIFNgVsControl')
    data_il17 = pd.read_excel(degs1_path, sheet_name='DEAnalysisIL17VsControl')
    data_tnfa = pd.read_excel(degs1_path, sheet_name='DEAnalysisTNFaVsControl')
    data_il13 = pd.read_excel(degs2_path, sheet_name='IL13')

    # Filter based on adj p val < adjp
    data_ifng_adjp = data_ifng[data_ifng['adj.P.Val'] < adjp]
    data_il17_adjp = data_il17[data_il17['adj.P.Val'] < adjp]
    data_tnfa_adjp = data_tnfa[data_tnfa['adj.P.Val'] < adjp]
    data_il13_adjp = data_il13[data_il13['adj.P.Val'] < adjp]

    # Filter based on lfc
    data_ifng_adjp = data_ifng_adjp[(data_ifng_adjp['logFC'] > lfc) | (data_ifng_adjp['logFC'] < -lfc)]
    data_il17_adjp = data_il17_adjp[(data_il17_adjp['logFC'] > lfc) | (data_il17_adjp['logFC'] < -lfc)]
    data_tnfa_adjp = data_tnfa_adjp[(data_tnfa_adjp['logFC'] > lfc) | (data_tnfa_adjp['logFC'] < -lfc)]
    data_il13_adjp = data_il13_adjp[(data_il13_adjp['logFC'] > lfc) | (data_il13_adjp['logFC'] < -lfc)]

    #print(data_ifng_adjp)

    # Save lists of degs only (ensembl ids)
    data_ifng_adjp.to_csv(ifng_degs_out, sep='\t', index=False, header=True)
    data_il17_adjp.to_csv(il17_degs_out, sep='\t', index=False, header=True)
    data_tnfa_adjp.to_csv(tnfa_degs_out, sep='\t', index=False, header=True)
    data_il13_adjp.to_csv(il13_degs_out, sep='\t', index=False, header=True)

    frames = [data_ifng_adjp, data_il17_adjp, data_tnfa_adjp, data_il13_adjp]

    return(frames)


def load_omnipath(omnipath_path):
    """
    Load and preprocess the OmniPath signalling prior network.

    :param omnipath_path: Path to the omnipath table (txt file)

    :return: list of filtered dataframes of degs ('frames')
    """
    # Load unfiltered omnipath data
    omnip = pd.read_csv(omnipath_path, sep="\t")

    omnip['source_uniprotac'] = omnip['source_uniprotac'].str.upper()
    omnip['target_uniprotac'] = omnip['target_uniprotac'].str.upper()
    omnip.rename(columns={'source_uniprotac':'source_uniprot','target_uniprotac':'target_uniprot'}, inplace=True)
    #omnip_min = omnip[['source_uniprotac', 'target_uniprotac']]
    #print(omnip)

    return(omnip)


def load_expressed(expressed1_path, expressed2_path, exp_filt):
    """
    Import the files of expresssion data and create lists of expressed genes

    :param expressed1_path: Path to the excel file with expression data from IL17, IFNg and TNFa (from Nick and Polychronis)

    :param expressed2_path: Path to the excel file with expression data from IL13 (from Nick and Polychronis)

    :param exp_filt: Cut off value for expression of gene (using >). At least 2 of 3 replicates need this fpkm value.

    :return: List of lists. Each list is expressed genes in one cytokine condition. in correct order based on df_names.

    """
    # Load the file/s
    data_1 = pd.read_excel(expressed1_path, sheet_name='fpkm')
    data_2 = pd.read_excel(expressed2_path, sheet_name='fpkm')

    # Get right columns
    il17_exp = data_1[['geneID', 'IL17_1', 'IL17_2', 'IL17_3', 'IL17_4']]
    ifng_exp = data_1[['geneID', 'IFNg_1', 'IFNg_2', 'IFNg_3', 'IFNg_4']]
    tnfa_exp = data_1[['geneID', 'TNFa_1', 'TNFa_2', 'TNFa_3', 'TNFa_4']]
    il13_exp = data_2[['geneID', 'SH_IL13', 'MW_IL13', 'ZB_IL13']]
    all_exp = [ifng_exp,il17_exp,tnfa_exp,il13_exp]

    # List of list of exp genes
    exp_lists = list()

    # Apply filters for expression

    for df in all_exp:
        df = df.set_index('geneID')
        #  > exp_filt in 2 or more cols
        c = (df > exp_filt).sum(axis=1)
        c = c[c >= 2]
        exp_genes = list(c.index.values)

        # Convert list to pandas series
        exp_series = pd.Series(exp_genes).astype(str)
        exp_series = exp_series.to_frame(name=None)

        # Create list of list of expressed genes for each cytokine in correct order
        exp_lists.append(exp_series)
        #print(exp_lists)

    # Save the series's of expressed genes
    for idx, expressed in enumerate(exp_lists):
        filepath = f"{df_names[idx]}_expressed_genes.txt"
        exp_out = f"{os.path.join(exp_out_dir, filepath)}"
        # Save the list of expressed genes to file
        expressed.to_csv(exp_out, sep='\t', index=False)

    return exp_lists


def network_to_ensembl(network, idconversion, database_type):
    """
    Convert omnipath and tftg networks to Ensembl Ids from uniprot or gene symbol

    :param network: Network with symbol ids in columns 'TF' and 'target'.

    :param idconversion: dataframe with columns 'Gene stable ID' and 'Gene name'. Long format.

    :return: network with ensembl ids and gene symbols in another column.

    """

    if database_type == "omnipath":
        idconversion_u = idconversion.groupby('UniProtKB/Swiss-Prot ID')['Gene stable ID'].agg(
            lambda col: ','.join(col)).reset_index()

        # Join uniprot ids to the network. Remove rows where no match found.
        network_e = network.merge(idconversion_u, left_on='source_uniprot', right_on='UniProtKB/Swiss-Prot ID',
                                  how='inner')
        network_e.drop(['UniProtKB/Swiss-Prot ID'], axis=1, inplace=True)
        network_e = network_e.rename(columns={'Gene stable ID': 'source_ensembl_id'})
        network_e2 = network_e.merge(idconversion_u, left_on='target_uniprot', right_on='UniProtKB/Swiss-Prot ID',
                                     how='inner')
        network_e2.drop(['UniProtKB/Swiss-Prot ID'], axis=1, inplace=True)
        network_e2 = network_e2.rename(columns={'Gene stable ID': 'target_ensembl_id'})

        # Expand the multi ens ids into multiple rows
        network_e3 = network_e2['target_ensembl_id'].str.split(',', expand=True).stack().reset_index(level=0).set_index(
            'level_0').rename(columns={0: 'target_ensembl_id'}).join(network_e2.drop('target_ensembl_id', 1))
        network_e3 = network_e3.reset_index(level=0)
        network_e4 = network_e3['source_ensembl_id'].str.split(',', expand=True).stack().reset_index(level=0).set_index(
            'level_0').rename(columns={0: 'source_ensembl_id'}).join(network_e3.drop('source_ensembl_id', 1))
        network_e4.drop(['index'], axis=1, inplace=True)
    else:
        # Make the id conversion table uniprot unique
        idconversion_u = idconversion.groupby('Gene name')['Gene stable ID'].agg(
            lambda col: ','.join(col)).reset_index()

        # Join uniprot ids to the network. Remove rows where no match found.
        network_e = network.merge(idconversion_u, left_on='TF', right_on='Gene name', how='inner')
        network_e.drop(['Gene name'], axis=1, inplace=True)
        network_e = network_e.rename(columns={'Gene stable ID': 'source_ensembl_id'})
        network_e2 = network_e.merge(idconversion_u, left_on='target', right_on='Gene name', how='inner')
        network_e2.drop(['Gene name'], axis=1, inplace=True)
        network_e2 = network_e2.rename(columns={'Gene stable ID': 'target_ensembl_id'})

        # Expand the multi ens ids into multiple rows
        network_e3 = network_e2['target_ensembl_id'].str.split(',', expand=True).stack().reset_index(level=0).set_index(
            'level_0').rename(columns={0: 'target_ensembl_id'}).join(network_e2.drop('target_ensembl_id', 1))
        network_e3 = network_e3.reset_index(level=0)
        network_e4 = network_e3['source_ensembl_id'].str.split(',', expand=True).stack().reset_index(level=0).set_index(
            'level_0').rename(columns={0: 'source_ensembl_id'}).join(network_e3.drop('source_ensembl_id', 1))
        network_e4.drop(['index'], axis=1, inplace=True)

    # Remove duplicate interactions
    network_e4 = network_e4.drop_duplicates(['target_ensembl_id','source_ensembl_id'], keep = 'first')

    return(network_e4)


def tftg_filtering(degs, tfnet):
    """
    Integrate the TFTGs with DEGs using ID conversion

    :param degs: dataframe of degs with all columns

    :param tfnet: dataframe of all possible tf regulatory edges from filtered input file (just source and target columns)

    :return: TF-> deg network where all tfs are degs also (source and target columns for uniprot and ensembl).
            Ensembl not necessarily unique.
    """

    # Extract the interesting columns of the deg df
    degs = degs[['ENSEMBL ID', 'adj.P.Val', 'logFC', 'SYMBOL']]

    # Join the deg table to the semi filtered tftg table, keep only rows where deg is in the source column
    tfnet_filtered = tfnet.merge(degs, left_on='source_ensembl_id', right_on='ENSEMBL ID',
                                           how='inner')
    tfnet_filtered = tfnet_filtered.drop('ENSEMBL ID', axis=1)
    tfnet_filtered = tfnet_filtered.rename(columns={'SYMBOL': 'source_symbol', 'adj.P.Val': 'source_adj.P.Val',
                                                                                             'logFC': 'source_logFC'})

    # Join the deg table to the tftg table, keep only rows where deg is in the target column
    tfnet_filtered2 = tfnet_filtered.merge(degs, left_on= 'target_ensembl_id', right_on= 'ENSEMBL ID', how='inner')
    tfnet_filtered2 = tfnet_filtered2.drop('ENSEMBL ID', axis=1)
    tfnet_filtered2 = tfnet_filtered2.rename(columns={'SYMBOL': 'target_symbol', 'adj.P.Val': 'target_adj.P.Val',
                                                    'logFC': 'target_logFC'})

    #print(tfnet_filtered2.columns.values)
    #print(tfnet.shape)

    # Output the deg<-TF networks
    return(tfnet_filtered2)


def filter_regulators(tf_deg_net, tf_filt_num):
    """
    Filter the tfs in the tf-deg networks by selecting those which target at least x degs (x = tf_filt_num)

    :param tf_deg_net: TF-> deg network where all tfs are degs also (source and target columns for ensembl).
    :param tf_filt_num: Minimum number of target nodes each regulator must have

    :return: Filtered TF-> deg networks where all TF target many degs.
    """

    #print(list(tf_deg_net))

    # Check there are no repeated interactions
    tf_deg_net = tf_deg_net.drop_duplicates()

    # Group the df by regulator to get counts of each regualtor in the network
    tf_net_regs = tf_deg_net.groupby(['source_ensembl_id'], as_index=False).count()

    # Get the regulators with count more than specified (tf_filt_num)
    tf_net_regs = tf_net_regs[tf_net_regs.target_ensembl_id >= tf_filt_num]

    # Filter the network for these regulators only
    #tf_net_filt = tf_deg_net[tf_deg_net['source_ensembl_id'].isin(tf_net_regs.source_ensembl_id)]
    tf_net_filt = tf_deg_net.loc[tf_deg_net.source_ensembl_id.isin(tf_net_regs.source_ensembl_id)].copy()

    # Add columns to say if the node is a tf, target or tf+target
    tf_net_filt['source_transcription_factors'] = "TF"
    tf_net_filt['target_tf_target'] = "TF_target"

    # Check the regulators have lots of targets
    # print(tf_net_filt.groupby(['source_uniprot'], as_index=False).count())

    return tf_net_filt


def omnipath_filtering(expression, omnip):
    """
    Filter the omnipath interactions so all nodes are expressed - except the cytokines which are exempt.
    Note that there is a chance the middle of the signalling pathway contains these cytokines now that I've
    removed the expression filter on them.

    :param expression: List/series of expressed genes (ensembl ID)
    :param omnip: Omnipath network to filter with source nodes in col 'source_ensembl_id' and target nodes in column
                target_ensembl_id'

    :return: Filtered omnipath network where all nodes are degs (or cytokines).

    """

    # List of omnipath ensembl IDs
    cytokine_ens = pd.Series(["ENSG00000111537","ENSG00000127318","ENSG00000112115","ENSG00000228849",
                               "ENSG00000228321", "ENSG00000206439", "ENSG00000204490", "ENSG00000232810",
                               "ENSG00000230108", "ENSG00000223952","ENSG00000169194"])

    # Add cytokine ids to expressed list
    exp_list = expression[0]
    #print(exp_list)
    exp_list = exp_list.append(cytokine_ens)

    # Join the expressed table to the ominpath table, keep only rows where expressed gene is in the source column
    omnip_filtered = omnip[(omnip['source_ensembl_id'].isin(exp_list)) & (omnip['target_ensembl_id'].isin(exp_list))]

    #print(len(omnip))
    #print(len(omnip_filtered))

    return omnip_filtered


def path_finding(omnipath, tfs, cytokine):
    """
    Get the shortest signalling paths between a cytokine and the TFs. Called for each cytokine separately.

    :param omnipath: Omnipath network with source nodes in col 'source_ensembl_id' and target nodes in column
            target_ensembl_id'.
    :param tfs: Series of TFs (ensembl ids) for cytokine of interest (ensembl ids)
    :param cytokine: String of the cytokine of interest eg. il17.

    :return: Network of all shortest paths from cytokine to tfs. Source nodes in col 'source'
            and target nodes in column 'target'.
    """

    # Convert the omnipath pandas df into a directed networkx network
    omni_x = nx.from_pandas_edgelist(omnipath, 'source_ensembl_id', 'target_ensembl_id', create_using=nx.DiGraph())

    cytokine_dict = {
        "ifng": "ENSG00000111537",
        "il17": "ENSG00000112115",
        "tnfa": "ENSG00000228321,ENSG00000206439,ENSG00000204490,ENSG00000232810,ENSG00000230108,ENSG00000223952",
        "il13": "ENSG00000169194"}

    # Empty dataframe for the edges
    shortest_net_all = pd.DataFrame()

    #print(type(omni_x))

    if (cytokine=="tnfa"):
        # Account for having >1 id
        value1 = cytokine_dict[cytokine].split(",")
        for val in value1:
            for tf in tfs:
                # Get all the shortest paths between source and target
                try:
                    for p in nx.all_shortest_paths(omni_x, source=val, target=tf):

                        if len(p) <= 6:  # Don't include where paths are longer than 5 connections

                            # Convert path to edges and append to df
                            for idx, item in enumerate(p):

                                if (len(p) - 1) == idx:
                                    continue
                                if idx == 0:
                                    level = "cytokine-receptor"
                                elif idx == 1:
                                    level = "receptor-interm1"
                                elif idx == (len(p) - 2):
                                    level = "final_interm-tf"
                                else:
                                    level = "interm-interm"
                                next = idx + 1
                                edges = {'source': [item],
                                         'target': [p[next]],
                                         'level': [level]}
                                df_edges = pd.DataFrame(edges, columns=['source', 'target','level'])
                                # Append the edges to the others for the cytokine
                                shortest_net_all = shortest_net_all.append(df_edges)
                except nx.NetworkXNoPath:
                    print("no paths found for " + cytokine + ", source= " + val + ", target = " + tf)
                    continue
    else:
        for tf in tfs:
            # Get all the shortest paths between source and target

            try:
                paths = nx.all_shortest_paths(omni_x, source=cytokine_dict[cytokine], target=tf)
                #print("Paths found for " + cytokine + ", source= " + cytokine_dict[cytokine] + ", target = " + tf)

                for p in paths:
                    if len(p) <= 6:  # Don't include where paths are longer than 5 connections
                        # Convert path to edges and append to df
                        for idx, item in enumerate(p):
                            if (len(p) - 1) == idx:
                                continue
                            if idx == 0:
                                level = "cytokine-receptor"
                            elif idx == 1:
                                level = "receptor-interm1"
                            elif idx == (len(p) - 2):
                                level = "final_interm-tf"
                            else:
                                level = "interm-interm"
                            next = idx + 1

                            edges = {'source': [item],
                                     'target': [p[next]],
                                     'level': [level]}
                            df_edges = pd.DataFrame(edges, columns=['source', 'target', 'level'])
                            # Append the edges to the others for the cytokine
                            shortest_net_all = shortest_net_all.append(df_edges)

            except nx.NetworkXNoPath:
                print("no paths found for " + cytokine + ", source= " + cytokine_dict[cytokine] + ", target = " + tf)
                continue

    if not shortest_net_all.empty:

        # Remove duplicates in shortest_net_all, joining the level info together
        shortest_net_all = shortest_net_all.groupby(['source', 'target'], as_index=False).agg(
            {'level': ','.join}).reset_index()

        # Append the other omnipath data columns
        #shortest_net_all2 = omnipath[(omnipath['source_ensembl_id'].isin(shortest_net_all["source"])) &
        #                              (omnipath['target_ensembl_id'].isin(shortest_net_all["target"]))]
        shortest_net_all2 = omnipath.merge(shortest_net_all, left_on=['source_ensembl_id','target_ensembl_id'],
                                right_on=['source','target'], how='inner')
        shortest_net_all2.drop(['source','target','index'], axis=1, inplace=True)
        #print(shortest_net_all2)
        return shortest_net_all2

    else:
        return shortest_net_all


def filter_hpa(hpa_path, network):
    """
    Filter network so that all nodes are expressed in the colon (human protein atlas)

    :param hpa_path: Path to the hpa colon dataset downloaded from hpa website
    :param network: Network to filter. Cols: source_ensembl_id, target_ensembl_id (and others)

    :return: Filtered network with same columns.

    """

    # Open hpa dataset
    hpa_df = pd.read_csv(hpa_path, sep = "\t")
    hpa_list = hpa_df['Ensembl']

    # Add cytokines to hpa list so we don't exclude them
    cytos = pd.Series(["ENSG00000111537","ENSG00000127318","ENSG00000112115","ENSG00000228849",
                               "ENSG00000228321", "ENSG00000206439", "ENSG00000204490", "ENSG00000232810",
                               "ENSG00000230108", "ENSG00000223952","ENSG00000169194"])
    hpa_list = hpa_list.append(cytos)

    # Filter using the ensembl id columns of the network
    hpa_network = network[network.source_ensembl_id.isin(hpa_list)]
    hpa_network = hpa_network[hpa_network.target_ensembl_id.isin(hpa_list)]

    #print(len(hpa_network))
    #print(len(network))

    return hpa_network


def tidy_omnipath(omnipath,degs):
    """
    Sort out the columns of the omnipath network so it matches the tf network

    :param omnipath: Omnipath network
    :param degs: The table of diff exp data with all columns

    :return: Omnipath network with tidied columns and deg data appended where present

    """


    omnipath['interaction_evidence_psimi'] = omnipath["psi-mi"].map(str) + ";" + omnipath["psi-mi2"].map(str) + ";" + omnipath[
        "psi-mi3"].map(str) + ";" + omnipath["psi-mi4"].map(str)
    omnipath.rename(columns={"pubmed": "interaction_evidence_pubmed", "taxid": "source_taxid"}, inplace=True)
    omnipath['target_taxid'] = omnipath['source_taxid']
    omnipath.drop(['psi-mi', 'psi-mi2', 'psi-mi3', 'psi-mi4'], axis=1, inplace=True)

    # Add diff exp data to the omnipath dataset (so it matches the TF dataset)
    omnipath1 = omnipath.merge(degs, left_on='source_ensembl_id', right_on='ENSEMBL ID',
                                 how='left')
    omnipath1 = omnipath1.drop(['ENSEMBL ID', 'zscore', 'P.Value', 'coef Control'], axis=1)
    omnipath1 = omnipath1.loc[:, ~omnipath1.columns.str.startswith('coef')]
    omnipath1 = omnipath1.rename(columns={'SYMBOL': 'source_symbol', 'adj.P.Val': 'source_adj.P.Val',
                                                    'logFC': 'source_logFC'})

    # Join the deg table to the tftg table, keep only rows where deg is in the target column
    omnipath2 = omnipath1.merge(degs, left_on='target_ensembl_id', right_on='ENSEMBL ID', how='left')
    omnipath2 = omnipath2.drop(['ENSEMBL ID','zscore','P.Value','coef Control'], axis=1)
    omnipath2 = omnipath2.loc[:, ~omnipath2.columns.str.startswith('coef')]
    omnipath2 = omnipath2.rename(columns={'SYMBOL': 'target_symbol', 'adj.P.Val': 'target_adj.P.Val',
                                                      'logFC': 'target_logFC'})
    omnipath2 = omnipath2

    #print(list(omnipath2))

    return (omnipath2)


######################
# Run the main script
######################

if __name__ == "__main__":

    # Call functions for loading the data (except omnipath - that is loaded later)
    # Convert tftg network to ensembl ids
    if do_all_load:
        idconversion, idconversion_sym = load_id_conversion(idconversion_path, idconversion_path_sym)
        tftgs = load_tftg_prior_knowledge(tftg_path)
        all_degs = load_degs(degs1_path, degs2_path)
        exp_lists = load_expressed(expressed1_path, expressed2_path, exp_filt)
        # Convert the tftg network to ensembl
        tftg_u = network_to_ensembl(tftgs, idconversion_sym, "tftgs")
        # and save
        tftg_u.to_csv(tftgs_out, sep='\t', index=False, header=True)
    else:
        # To speed subsequent running load the data which was saved during the load functions.

        # Load id conversion table
        idconversion = pd.read_csv(idconversion_out, sep='\t')
        idconversion_sym = pd.read_csv(idconversion_sym_out, sep='\t')

        # Load tftgs table (with ensembl ids)
        tftg_u = pd.read_csv(tftgs_out, sep='\t')

        # Load lists of degs
        ens_data_ifng_adjp = pd.read_csv(ifng_degs_out, sep='\t')
        ens_data_il17_adjp = pd.read_csv(il17_degs_out, sep='\t')
        ens_data_tnfa_adjp = pd.read_csv(tnfa_degs_out, sep='\t')
        ens_data_il13_adjp = pd.read_csv(il13_degs_out, sep='\t')
        all_degs =  [ens_data_ifng_adjp, ens_data_il17_adjp,ens_data_tnfa_adjp, ens_data_il13_adjp]

        # Load the lists of expressed genes
        exp_lists = []
        for cytok in df_names:
            filename = f"{cytok}_expressed_genes.txt"
            filepath = f"{os.path.join(exp_out_dir, filename)}"
            # Open the list of expressed genes file
            exp_list = pd.read_csv(filepath, sep='\t', header=None)
            exp_lists.append(exp_list)

    # Call functions for generating tf-deg networks
    tf_deg_nets = []
    if do_tfnet:
       # Go through each deg list, get tf-deg networks (filter so all regs are degs also) and save to file
        for idx, deg_list in enumerate(all_degs):

            # Run function to get the TF-> DEG networks and filter so all nodes are degs
            tfnet_filtered1 = tftg_filtering(deg_list, tftg_u)
            # If required, filter the network to ensure all nodes are expressed in the colon (human protein atlas)
            if do_filter_hpa:
                tfnet_filtered2 = filter_hpa(hpa_path, tfnet_filtered1)
            else:
                tfnet_filtered2 = tfnet_filtered1
            # Specify output directory
            filepath = f"{df_names[idx]}_tf_deg_network.txt"
            deg_tf_out = f"{os.path.join(deg_tf_output_dir, filepath)}"
            # Save the deg<-TF network to file
            tfnet_filtered2.to_csv(deg_tf_out, sep='\t', index=False, header=True)
            # Append the network to list of networks
            tf_deg_nets.append(tfnet_filtered2)
    else:
        # Load pre-saved tf-deg networks
        for name in df_names:
            # Get filepaths of the networks
            filepath1 = f"{name}_tf_deg_network.txt"
            filepath2 = f"{os.path.join(deg_tf_output_dir, filepath1)}"
            tfnet_filtered1 = pd.read_csv(filepath2, sep='\t')
            tf_deg_nets.append(tfnet_filtered1)

    # Call function to filter tf-deg networks for important regulators using degree
    tf_deg_nets2 = []
    if do_tf_filt:
        # Filter the tfs in the tf-deg networks by selecting those which target at least x degs (x = tf_filt_num)
        for idx, tf_net in enumerate(tf_deg_nets):
            # Run function to do the filter
            tf_net_filt1 = filter_regulators(tf_net, tf_filt_num)
             # Specify output directory
            filepath = f"{df_names[idx]}_tf_deg_network_regfilt.txt"
            deg_tf_filt_out = f"{os.path.join(deg_tf_filt_output_dir, filepath)}"
            # Save the deg<-TF regulator filtered network to file
            tf_net_filt1.to_csv(deg_tf_filt_out, sep='\t', index=False, header=True)

            # Get list of regulator filtered tf-deg dataframes
            tf_deg_nets2.append(tf_net_filt1)
        tf_deg_nets = tf_deg_nets2

    # Call function to filter tf-deg networks for important regulators using degree and lfc (CHAT)
    tf_deg_nets2 = []
    if do_tf_filt_chat:
        # Otherwise filter regs using chat (outside of this script) and then load in the filtered networks
        for name in df_names:
            # Specify file path
            filepath = f"{name}_tf_deg_network_chatfilt_{chat_pcorr}.txt"
            filepath2 = f"{os.path.join(chat_path, filepath)}"
            # Open the file
            tf_net_filt1 = pd.read_csv(filepath2, sep='\t')
            # Append to list of dataframes
            tf_deg_nets2.append(tf_net_filt1)
        tf_deg_nets = tf_deg_nets2

    # Start omnipath data analysis and joining of omipath to tf networks
    if do_omnipath:
        # Load the omnipath data
        omnip1 = load_omnipath(omnipath_path)

        # Convert omnipath to ensembl IDs
        omnip2 = network_to_ensembl(omnip1,idconversion,"omnipath")

        # Empty list for each cytokine network to be listed
        all_cyto_nets = []
        # Iterate cytokines
        for idx, deg_list in enumerate(all_degs):

            omnipath0 = omnip2.copy()
            # Sort out the columns of omnipath data
            omnip3 = tidy_omnipath(omnipath0, deg_list)

            # Carry out expression filter of omnipath genes (except the cytokines themselves)
            if do_omnipath_filter_exp:

                # Filter omnipath for expressed genes/proteins only
                omnip_filtered1 = omnipath_filtering(exp_lists[idx], omnip3)

                # Save the deg filtered omnipath networks
                #filepath = f"{df_names[idx]}_omnipath_expressed_network.txt"
                #deg_omnip_out = f"{os.path.join(exp_omnip_out_dir, filepath)}"
                #omnip_filtered1.to_csv(deg_omnip_out, sep='\t', index=False, header=True)

            else:
                omnip_filtered1 = omnip3

            # Carry out the HPA-colon filter of omnipath genes (except the cytokines themselves)
            if do_filter_hpa:
                omnip_filtered2 = filter_hpa(hpa_path, omnip_filtered1)

                # Save the hpa filtered omnipath networks
                #filepath = f"{df_names[idx]}_omnipath_expressed_network.txt"
                #deg_omnip_out = f"{os.path.join(exp_omnip_out_dir, filepath)}"
                #omnip_filtered1.to_csv(deg_omnip_out, sep='\t', index=False, header=True)
            else:
                omnip_filtered2 = omnip_filtered1

            # Get the list of TFs
            tf_net = tf_deg_nets[idx]
            uniq_tfs = pd.Series(tf_net['source_ensembl_id']).unique()

            # Get the shortest path between the cytokine and the TFs (in list)
            shortest_nets_all1 = path_finding(omnip_filtered2, uniq_tfs, df_names[idx])

            if not shortest_nets_all1.empty:
                # Append cytokine->tf network to list
                # all_cyto_nets.append(shortest_nets_all1)

                # Save the cytokine->tf networks to file
                filepath = f"{df_names[idx]}_omnip_tf_network.txt"
                filepath2 = f"{os.path.join(omnip_tf_out_dir, filepath)}"
                shortest_nets_all1.to_csv(filepath2, sep='\t', index=False, header=True)

                # Get list of tfs targetted by the cytokines
                tflevel = shortest_nets_all1[shortest_nets_all1.level.str.contains("final_interm-tf")]
                targeted_tfs = pd.Series(tflevel['target_ensembl_id']).unique()

                # Filter the tf network so all tfs are targetted by cytokine
                tf_net2 = tf_net[tf_net.source_ensembl_id.isin(targeted_tfs)]

                # Join the tf network to the cytokine/omnipath network
                tf_net2['level'] = "tf-target"
                all_net = tf_net2.append(shortest_nets_all1)
                all_net = all_net.drop_duplicates()

                # Save output
                filepath = f"{df_names[idx]}_omnip_deg_network.txt"
                filepath2 = f"{os.path.join(omnip_deg_out_dir, filepath)}"
                all_net.to_csv(filepath2, sep='\t', index=False, header=True)
