# First, ensure you have the necessary libraries installed:
# pip install rdkit-pypi networkx

import sys
try:
    from rdkit import Chem
    import networkx as nx
except ImportError:
    print("Error: RDKit and NetworkX libraries are required.")
    print("Please install them by running: pip install rdkit-pypi networkx")
    sys.exit(1)

# Step 1: Identify the molecule and represent it.
# The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
# We will use its SMILES string and include all hydrogen atoms as requested.
smiles = 'c1cc2c3ccc4cccc5c4c(c2c(c1)S)c3cc5'
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# Step 2: Convert the molecule into a graph representation.
G = nx.Graph()
for atom in mol.GetAtoms():
    G.add_node(atom.GetIdx())
for bond in mol.GetBonds():
    G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

# Step 3: Calculate the shortest path distances between all pairs of atoms.
path_lengths = dict(nx.all_pairs_shortest_path_length(G))

# Step 4: Calculate the Wiener Index (W).
# W is the sum of distances between all unique pairs of nodes.
wiener_index = 0
nodes = list(G.nodes())
for i in range(len(nodes)):
    for j in range(i + 1, len(nodes)):
        wiener_index += path_lengths[nodes[i]][nodes[j]]
        
# Step 5: Calculate the Szeged Index (Sz).
# For each edge (u, v), count the number of nodes closer to u (n_u) and closer to v (n_v).
# Sz is the sum of n_u * n_v over all edges.
szeged_index = 0
for u, v in G.edges():
    n_u = 0
    n_v = 0
    for w in G.nodes():
        dist_u = path_lengths[w][u]
        dist_v = path_lengths[w][v]
        if dist_u < dist_v:
            n_u += 1
        elif dist_v < dist_u:
            n_v += 1
    szeged_index += n_u * n_v

# Step 6: Calculate the ratio and print the final equation.
if wiener_index > 0:
    ratio = szeged_index / wiener_index
    # The final output is an equation showing the numbers involved.
    print(f"For the major reduction product, perylene-3-thiol (including H):")
    print(f"Szeged Index (Sz) = {szeged_index}")
    print(f"Wiener Index (W) = {wiener_index}")
    print(f"Final Equation for the Sz/W Ratio:")
    print(f"{szeged_index} / {wiener_index} = {ratio}")
else:
    print("Could not calculate the ratio as the Wiener index is zero.")
