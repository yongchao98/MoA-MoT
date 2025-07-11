import rdkit
from rdkit import Chem
import networkx as nx

def calculate_indices():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1: Define the molecule (perylene-3-thiol) using its SMILES string.
    # This is the major reduction product of di(perylene-3-yl) disulfide.
    smiles = "c1cc(S)c2c3c(c1)c1cccc4ccc5cccc(c14)c5c23"
    
    # Step 2: Create an RDKit molecule object and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    mol_with_hs = Chem.AddHs(mol)

    # Step 3: Convert the RDKit molecule to a NetworkX graph.
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    num_atoms = G.number_of_nodes()
    
    # Pre-calculate all-pairs shortest path lengths, which is needed for both indices.
    # This returns a dictionary of dictionaries: {source: {target: length}}.
    all_pairs_paths = dict(nx.all_pairs_shortest_path_length(G))

    # Step 4: Calculate the Wiener Index (W).
    # W = 0.5 * sum of all shortest path distances.
    wiener_index = 0
    for source in all_pairs_paths:
        for target in all_pairs_paths[source]:
            wiener_index += all_pairs_paths[source][target]
    wiener_index //= 2  # Each pair is counted twice, so divide by 2.

    # Step 5: Calculate the Szeged Index (Sz).
    # Sz = sum over all edges (u,v) of n_u(e) * n_v(e).
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0  # Nodes closer to u than to v
        n_v = 0  # Nodes closer to v than to u
        for w in G.nodes():
            dist_u_w = all_pairs_paths[u][w]
            dist_v_w = all_pairs_paths[v][w]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v
        
    # Step 6: Compute the ratio and print the results.
    ratio = szeged_index / wiener_index if wiener_index != 0 else 0

    print("Target Molecule: Perylene-3-thiol (C20H12S)")
    print(f"Total atoms (nodes in graph): {num_atoms}")
    print("\n--- Topological Index Calculation ---")
    print(f"Szeged Index (Sz): {szeged_index}")
    print(f"Wiener Index (W): {wiener_index}")
    print("\n--- Szeged/Wiener Index Ratio ---")
    print(f"{szeged_index} / {wiener_index} = {ratio}")
    
    # Final answer in the specified format
    # Using round() for a clean final answer.
    print(f"\n<<<{round(ratio, 4)}>>>")

if __name__ == "__main__":
    calculate_indices()