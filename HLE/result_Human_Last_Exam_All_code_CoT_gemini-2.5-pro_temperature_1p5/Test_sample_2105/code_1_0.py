import networkx as nx
from rdkit import Chem

def calculate_indices():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # 1. Define the molecule: perylene-3-thiol, the reduction product.
    # SMILES string for perylene-3-thiol
    smiles = 'Sc1cc2cccc3c2c4c(c1)ccc5cccc4c35'
    mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogen atoms as specified
    mol_with_hs = Chem.AddHs(mol)

    # 2. Convert the RDKit molecule to a NetworkX graph
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    num_atoms = mol_with_hs.GetNumAtoms()
    print(f"The molecule under analysis is perylene-3-thiol (C20H12S).")
    print(f"The molecular graph includes {num_atoms} atoms (vertices), including hydrogens.")
    
    # Pre-calculate all-pairs shortest paths, which is needed for both indices
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # 3. Calculate Wiener Index (W)
    wiener_index = 0
    # Sum of all distances between unordered pairs
    for source_node in path_lengths:
        for target_node in path_lengths[source_node]:
             # Add path length, divide total by 2 later to correct for double counting
            wiener_index += path_lengths[source_node][target_node]
    wiener_index //= 2
    
    # 4. Calculate Szeged Index (Sz)
    szeged_index = 0
    for u, v in G.edges():
        # For each edge (u, v), find the number of vertices closer to u or v
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

    # 5. Compute and print the ratio
    if wiener_index == 0:
        ratio = 0
    else:
        ratio = szeged_index / wiener_index

    print("\n--- Calculation Results ---")
    print(f"Szeged Index (Sz): {szeged_index}")
    print(f"Wiener Index (W): {wiener_index}")
    print(f"Sz/W Ratio = {szeged_index} / {wiener_index} = {ratio}")

if __name__ == '__main__':
    calculate_indices()