# First, we generate the Protein Structure Network based on the
# distances between the centers of mass of the side chains.
# All residues apart from glycine (that has no side chain) are used.
# A distance cut-off of 5.125 A is used.

sh 0-contact.sh > log_0

# Then, we filter the network and we calculate hubs and connected
# components.
# A minimum edge persistence value of 20.0 is used for filtering,
# meaning that only edges present in >= 20% of the frames of the
# simulation are kept in the filtered network.

sh 1-graph_analysis.sh > log_1

# Finally, we calculate the shortest paths between selected pairs
# of residues. Such calculations are performed inside each subfolder
# of the "paths" folder, by running (inside the subfolder):

sh get_paths.sh > log_paths

# Pairs of residues considered for shortest paths' calculation are:
# Y22-L108, Y22-H113, Y22-F146
# S23-L108, S23-H113, S23-F146
# Q26-L108, Q26-H113, Q26-F146
# S28-L108, S28-H113, S28-F146
# V155-L108, V155-H113, V155-F146

# Since the PDB numbering of the reference structure differs from
# the UniProt numbering because of the lack of the flexible loop
# between helices alpha1 and alpha2 in Bcl-xL (residues 45-84),
# residues starting from position 85 are numbered differently in
# paths. Specifically, you should add +40 to any residue number
# after 44 to retrieve the UniProt number of that residue.
