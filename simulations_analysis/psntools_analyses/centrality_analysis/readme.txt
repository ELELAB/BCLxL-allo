# In order to run the centrality analysis, run:

python3.7 centrality_analysis.py

# Reference structures are taken from: ../../reference_structures_psntools

# Files and folders needed to run:
# * node_labels.txt -> text file where the labels to be used for the nodes of the PSNs are stored.
# * selected_nodes.txt -> text file containing the nodes that should be included in the plots used
#                         in Figure 5 of the paper.
# * config_plot_files -> folder containing YAML configuration files for the plots to be generated,
#                        namely the heatmaps showing the betweenness centrality, closeness centrality
#                        or degree of each node of the PSNs/of selected nodes of the PSNs.

# The run will produce:
# * CSVs -> a folder containing the CSV files with the centrality values for all nodes of all PSNs.
# * plots_analysis -> a folder containing the PDF files of the heatmaps of the centrality values
#		      for all nodes of the PSNs.
# * plots_paper -> a folder containing the PDF files of the heatmaps of the centrality values for
#                  selected nodes of the PSNs. Those are the heatmaps shown in Figure 5 of the paper.
