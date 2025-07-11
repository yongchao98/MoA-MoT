# First, install the necessary libraries if you haven't already:
# pip install rdkit-pypi networkx

import networkx as nx
from rdkit import Chem

def solve_molecular_ratio():
    """
    This function calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1: Create a perylene molecule and add hydrogens using RDKit.
    # The SMILES string represents the perylene carbon skeleton.
    perylene_smiles = 'c1ccc2c3c(c1)c4ccc5ccccc5c4c3c2'
    mol = Chem.MolFromSmiles(perylene_smiles)
    mol_h = Chem.AddHs(mol)

    # Step 2: Convert the RDKit molecule object into a NetworkX graph.
    # Atoms become nodes, and bonds become edges.
    G = nx.Graph()
    for atom in mol_h.GetAtoms():
        G.add_node(atom.GetIdx(), symbol=atom.GetSymbol())

    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Identify the IUPAC C3 position (a "bay" carbon) and its attached hydrogen.
    # Bay carbons in perylene are bonded to 3 other carbons and 1 hydrogen.
    target_c_idx = -1
    old_h_idx = -1
    for atom in mol_h.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
            if len(c_neighbors) == 3 and len(h_neighbors) == 1:
                target_c_idx = atom.GetIdx()
                old_h_idx = h_neighbors[0].GetIdx()
                break # Found one bay carbon; they are all equivalent by symmetry.

    if target_c_idx == -1:
        print("Error: Could not find a target carbon (bay position) to modify.")
        return

    # Step 4: Modify the graph to form perylene-3-thiol.
    # Remove the original hydrogen atom.
    G.remove_node(old_h_idx)

    # Add new nodes for the thiol group (-SH).
    # We use high index numbers to avoid collisions.
    s_idx = mol_h.GetNumAtoms()
    h_new_idx = mol_h.GetNumAtoms() + 1
    
    G.add_node(s_idx, symbol='S')      # Add Sulfur atom
    G.add_node(h_new_idx, symbol='H') # Add Hydrogen atom for the thiol

    # Add new bonds for the thiol group (C-S-H).
    G.add_edge(target_c_idx, s_idx)
    G.add_edge(s_idx, h_new_idx)

    # Step 5: Calculate the Wiener and Szeged indices.
    # First, compute all-pairs shortest paths, which is needed for both indices.
    paths = dict(nx.all_pairs_shortest_path_length(G))

    # Calculate Wiener Index (W)
    wiener_index = 0
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            wiener_index += paths[nodes[i]][nodes[j]]
    
    # Calculate Szeged Index (Sz)
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in G.nodes():
            dist_u_w = paths[w][u]
            dist_v_w = paths[w][v]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 6: Compute and print the final ratio.
    ratio = szeged_index / wiener_index
    
    print(f"Molecule: Perylene-3-thiol")
    print(f"Total number of atoms (nodes): {G.number_of_nodes()}")
    print("-" * 30)
    print(f"Wiener Index (W): {int(wiener_index)}")
    print(f"Szeged Index (Sz): {szeged_index}")
    print("-" * 30)
    print(f"Szeged/Wiener Index Ratio Calculation:")
    print(f"{szeged_index} / {int(wiener_index)} = {ratio}")

if __name__ == '__main__':
    solve_molecular_ratio()