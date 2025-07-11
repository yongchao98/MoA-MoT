import networkx as nx
try:
    from rdkit import Chem
except ImportError:
    print("RDKit not found. Please install it using: pip install rdkit-pypi")
    exit()

def calculate_indices():
    """
    Calculates the Wiener and Szeged indices for perylene-3-thiol and their ratio.
    """
    # Step 1: Define the molecule (perylene-3-thiol) from its SMILES string.
    # The SMILES for perylene-3-thiol is C1=CC=C2C3=C4C(C=CC5=C4C(C2=C1)=CC=C5)=C(C=C3)S.
    # The question specifies "including H", so we add explicit hydrogens.
    smiles = 'C1=CC=C2C3=C4C(C=CC5=C4C(C2=C1)=CC=C5)=C(C=C3)S'
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print("Error: Could not create molecule from the provided SMILES string.")
        return

    mol_h = Chem.AddHs(mol)

    # Step 2: Convert the RDKit molecule to a NetworkX graph.
    G = nx.Graph()
    for atom in mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Calculate all-pairs shortest paths, which are needed for both indices.
    # This returns a dictionary where path_lengths[u][v] is the shortest distance.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # Step 4: Calculate the Wiener Index (W).
    # It's the sum of shortest distances over all unique pairs of nodes.
    wiener_index = 0
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            wiener_index += path_lengths[nodes[i]][nodes[j]]

    # Step 5: Calculate the Szeged Index (Sz).
    # Sz = sum over all edges {e=(u,v)} of (n_u * n_v)
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in G.nodes():
            dist_u_w = path_lengths[u][w]
            dist_v_w = path_lengths[v][w]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 6: Compute the ratio and print the final results.
    if wiener_index > 0:
        ratio = szeged_index / wiener_index
        print(f"The analysis is for perylene-3-thiol (C20H12S), the reduction product.")
        print(f"Total number of atoms (graph nodes), including hydrogen: {G.number_of_nodes()}")
        print("-" * 40)
        print(f"Calculated Wiener Index (W) = {wiener_index}")
        print(f"Calculated Szeged Index (Sz) = {szeged_index}")
        print("-" * 40)
        # Final output showing the equation as requested
        print(f"Ratio (Sz / W) = {szeged_index} / {wiener_index} = {ratio}")
    else:
        print("Calculation failed: The Wiener index is zero or invalid.")

if __name__ == "__main__":
    calculate_indices()