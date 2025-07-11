# First, ensure you have the necessary libraries installed:
# pip install rdkit networkx numpy

import networkx as nx
from rdkit import Chem

def solve_molecular_indices():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1 & 2: Define the molecule (perylene-3-thiol) and create its graph representation
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    # We use its SMILES string to generate the molecule.
    smiles = "c1cc2c3c(c(S)cc4c3ccc5cccc1c54)cc2"
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogens to the molecule as specified in the problem
    mol = Chem.AddHs(mol)

    if mol is None:
        print("Error: Could not create molecule from SMILES string.")
        return

    # Create a NetworkX graph from the molecule
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    graph = nx.from_numpy_array(adj_matrix)

    # Step 3: Calculate all-pairs shortest paths
    # This is a dictionary of dictionaries: path_lengths[source][target] = length
    path_lengths = dict(nx.all_pairs_shortest_path_length(graph))

    # Step 4: Calculate the Wiener Index (W)
    wiener_index = 0
    nodes = list(graph.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            # Sum of shortest paths between all unique pairs of nodes
            wiener_index += path_lengths[nodes[i]][nodes[j]]

    # Step 5: Calculate the Szeged Index (Sz)
    szeged_index = 0
    for u, v in graph.edges():
        n_u = 0
        n_v = 0
        for w in graph.nodes():
            dist_wu = path_lengths[w][u]
            dist_wv = path_lengths[w][v]
            if dist_wu < dist_wv:
                n_u += 1
            elif dist_wv < dist_wu:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 6: Compute the ratio
    if wiener_index == 0:
        ratio = 0.0
    else:
        ratio = szeged_index / wiener_index

    # Print the results in a clear format
    print("Molecule: Perylene-3-thiol (including H)")
    print(f"Number of atoms (vertices): {graph.number_of_nodes()}")
    print(f"Number of bonds (edges): {graph.number_of_edges()}")
    print("-" * 30)
    print(f"The Wiener index (W) is: {wiener_index}")
    print(f"The Szeged index (Sz) is: {szeged_index}")
    print("-" * 30)
    print(f"The Szeged/Wiener index ratio is: {szeged_index} / {wiener_index} = {ratio}")

    # Return the final answer in the requested format
    print(f"\n<<<{ratio}>>>")

if __name__ == "__main__":
    solve_molecular_indices()
