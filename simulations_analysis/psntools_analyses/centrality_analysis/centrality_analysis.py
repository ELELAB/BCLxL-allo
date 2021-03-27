#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-


# standard library
import os
# third-party packages
import MDAnalysis as mda
from psntools.PSN import PSN, PSNGroup
import psntools.analysis as psnanalysis
import psntools.plotting as psnplotting



# path to where PSN matrices are stored 
matrices_path = "../../../{:s}/pyinteraph/contact/5.125/hc-graph_filt.dat"
# path to where PDB files are stored
pdbs_path = "../../../reference_structures_psntools/{:s}.pdb"

# file containing the nodes of interest
f_selnodes = "selected_nodes.txt"
# nodes of interest
selected_nodes = [l.rstrip("\n") for l in open(f_selnodes, "r")]

# file containing the labels for the nodes of interest
f_nodelabels = "node_labels.txt"
# nodes of interest
node_labels = [l.rstrip("\n") for l in open(f_nodelabels, "r")]

# current working directory
cwd = os.getcwd()
# output directory where to store CSV files
out_csv_dir = os.path.join(cwd, "CSVs")
os.makedirs(out_csv_dir, exist_ok = True)
# output directory where to store plots for analysis
out_plot_analysis_dir = os.path.join(cwd, "plots_analysis")
os.makedirs(out_plot_analysis_dir, exist_ok = True)
# output directory where to store plots for the paper
out_plot_paper_dir = os.path.join(cwd, "plots_paper")
os.makedirs(out_plot_paper_dir, exist_ok = True)

# directory containing customized configuration 
# YAML files for plotting
config_plot_dir = os.path.join(cwd, "config_plot_files")

# names of the files without extension
filenames = ["free_2nm_tip3p", \
             "free_cap_tip3p", \
             "free_2nm_tips3p", \
             "bclxl_puma_3.5nm", \
             "bclxl_puma_capped"]

# labels to be used for the PSNs
labels = [r"Bcl-$\mathregular{x_L}$ (1)", \
          r"Bcl-$\mathregular{x_L}$ (2)", \
          r"Bcl-$\mathregular{x_L}$ (3)", \
          r"Bcl-$\mathregular{x_L}$-PUMA (1)", \
          r"Bcl-$\mathregular{x_L}$-PUMA (2)"]

# PSN matrices
matrices = [matrices_path.format(fn) for fn in filenames]

# pdb files of the systems the PSNs were calculated on
pdbs = [pdbs_path.format(fn) for fn in filenames]

# Universe objects
universes = [mda.Universe(pdb, pdb) for pdb in pdbs]

# centrality metrics to be calculated
metrics = [{"degree" : {}}, \
           {"betweenness_centrality" : {}}, \
           {"closeness_centrality" : {"wf_improved" : True}}]

# create the PSNgroup
group = PSNGroup(matrices = matrices, \
                 universes = universes, \
                 labels = filenames)

# for each metric of interest
for metric in metrics:

    # get the metric name
    metric_name = list(metric.keys())[0]

    # get the dataframe containing the value of the metric
    # for all nodes of all PSNs in the group
    df = psnanalysis.get_nodes_df_psngroup(psngroup = group, \
                                           metric = metric)

    # save the dataframe in a CSV file
    df.to_csv(os.path.join(out_csv_dir, metric_name + ".csv"), \
              sep = ",")
    
    # plot a heatmap for visual analysis with all nodes and
    # annotations (use the default YAML file, we do not
    # require specific aesthetics)
    outfile_analysis = os.path.join(out_plot_analysis_dir, \
                                    metric_name + ".pdf")
    configfile_analysis = f"heatmap_nodes_{metric_name}"
    psnplotting.plot_heatmap_nodes(\
                        df = df, \
                        outfile = outfile_analysis, \
                        configfile = configfile_analysis, \
                        nodes_per_page = 20, \
                        psn_labels = labels, \
                        node_labels = node_labels)

    # plot e heatmap to be used as a panel in figure 5 in the paper
    # (use a customized YAML file with modified aesthetics)
    outfile_paper = os.path.join(out_plot_paper_dir, \
                                 "hubs_" + metric_name + ".pdf")
    configfile_paper = os.path.join(config_plot_dir, \
                                    f"heatmap_nodes_{metric_name}.yaml")
    psnplotting.plot_heatmap_nodes(\
            df = df, \
            outfile = outfile_paper, \
            configfile = configfile_paper, \
            selected_nodes = selected_nodes, \
            nodes_per_page = len(selected_nodes), \
            psn_labels = labels, \
            node_labels = node_labels)

