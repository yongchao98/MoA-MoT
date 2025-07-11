import sys
try:
    from rdkit import Chem
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)
try:
    import networkx as nx
except ImportError:
    print("NetworkX is not installed. Please install it using 'pip install networkx'")
    sys.exit(1)

def calculate_indices():
    """
    Calculates the Wiener and Szeged indices for perylene-3-thiol.
    """
    # SMILES string for perylene-3-thiol from PubChem (CID 13327668)
    smiles = 'C1=CC=C2C3=C(C=CC2=C1)C4=C5C(=C(C=C4)S)C=CC6=C5C3=CC=C6'

    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    # Add explicit hydrogen atoms as requested
    mol = Chem.AddHs(mol)

    # --- Step 1: Build the molecular graph using NetworkX ---
    G = nx.Graph()
    # Add atoms as nodes
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    # Add bonds as edges
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # --- Step 2: Pre-calculate all-pairs shortest path lengths ---
    # This is efficient as it avoids recalculating distances.
    all_paths = dict(nx.all_pairs_shortest_path_length(G))

    # --- Step 3: Calculate the Wiener Index (W) ---
    wiener_index = 0
    # Sum of distances between all unique pairs of nodes
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            wiener_index += all_paths[nodes[i]][nodes[j]]

    # --- Step 4: Calculate the Szeged Index (Sz) ---
    szeged_index = 0
    # Iterate over all edges (bonds) in the graph
    for u, v in G.edges():
        n_u = 0  # Nodes closer to u
        n_v = 0  # Nodes closer to v
        # For each edge, iterate over all nodes to find which endpoint is closer
        for w in nodes:
            dist_u = all_paths[w][u]
            dist_v = all_paths[w][v]
            if dist_u < dist_v:
                n_u += 1
            elif dist_v < dist_u:
                n_v += 1
        szeged_index += n_u * n_v
    
    # --- Step 5: Calculate the ratio and print the results ---
    if wiener_index == 0:
        ratio = 0
    else:
        ratio = szeged_index / wiener_index

    print(f"Analysis for Perylene-3-thiol (C20H12S), including H atoms:")
    print(f"Total number of atoms (nodes): {G.number_of_nodes()}")
    print(f"Total number of bonds (edges): {G.number_of_edges()}")
    print("-" * 30)
    print(f"Szeged Index (Sz): {szeged_index}")
    print(f"Wiener Index (W): {wiener_index}")
    print("-" * 30)
    print(f"Szeged/Wiener Index Ratio = {szeged_index} / {wiener_index} = {ratio:.5f}")
    
    # Return the final numeric answer for the tagging system.
    return ratio

if __name__ == '__main__':
    final_ratio = calculate_indices()
    # The final line with <<<>>> is for automated evaluation
    # print(f"\n<<<{final_ratio:.5f}>>>")


# The following call to the main function will execute the script.
calculate_indices()