import networkx as nx
from rdkit import Chem

def calculate_molecular_indices():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol, including hydrogens.
    """
    # Step 1 & 2: Define the molecule and construct the graph.
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    # Its SMILES string is used to generate the molecular graph.
    smiles = "c1cc2c3c(c(S)c4ccccc4c3cc1)cc5ccccc25"
    
    # Create an RDKit molecule object from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Could not create molecule from the provided SMILES string.")
        return

    # Add explicit hydrogen atoms to the molecule as required.
    mol_h = Chem.AddHs(mol)
    
    # Create an empty NetworkX graph.
    G = nx.Graph()
    
    # Add atoms from the molecule as nodes in the graph.
    G.add_nodes_from(range(mol_h.GetNumAtoms()))
    
    # Add bonds from the molecule as edges in the graph.
    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    print(f"The analysis is for perylene-3-thiol (C20H12S).")
    print(f"The molecular graph contains {mol_h.GetNumAtoms()} atoms (vertices), including hydrogens.")
    
    # Step 3: Calculate Wiener Index (W).
    # This is the sum of shortest path distances between all pairs of nodes.
    wiener_index = nx.wiener_index(G)
    
    # Step 4: Calculate Szeged Index (Sz).
    szeged_index = 0
    # Pre-calculate all-pairs shortest path lengths for efficiency.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # Iterate over each edge (bond) in the graph.
    for u, v in G.edges():
        n_u = 0  # Counter for nodes closer to u
        n_v = 0  # Counter for nodes closer to v
        
        # For each node w, determine if it's closer to u or v.
        for w in G.nodes():
            dist_wu = path_lengths[w][u]
            dist_wv = path_lengths[w][v]
            if dist_wu < dist_wv:
                n_u += 1
            elif dist_wv < dist_wu:
                n_v += 1
        
        # Add the product of the counts to the Szeged index.
        szeged_index += n_u * n_v

    # Step 5: Compute and print the ratio.
    if wiener_index == 0:
        print("Error: Wiener index is zero, so the ratio cannot be calculated.")
        return

    ratio = szeged_index / wiener_index
    
    print("\n--- Calculation Result ---")
    print(f"Szeged index (Sz) / Wiener index (W) = {szeged_index} / {wiener_index} = {ratio}")

if __name__ == '__main__':
    try:
        calculate_molecular_indices()
    except ImportError:
        print("Required libraries are not installed. Please install them using:")
        print("pip install rdkit-pypi networkx")
