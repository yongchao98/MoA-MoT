import networkx as nx
from rdkit import Chem

def solve_szeged_wiener_ratio():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol including hydrogen atoms.
    """
    print("Step 1: Identifying the molecule...")
    print("The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.")
    
    # The canonical SMILES for perylene-3-thiol
    smiles = "c1cc(S)c2c3ccc4cccc5cccc(c1)c3c2c45"
    
    print("\nStep 2: Generating the molecular graph including all hydrogen atoms...")
    # Create RDKit molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles)
    # Add explicit hydrogen atoms
    mol_with_hs = Chem.AddHs(mol)
    
    # Convert RDKit molecule to NetworkX graph
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
    print(f"Graph created with {G.number_of_nodes()} nodes (atoms) and {G.number_of_edges()} edges (bonds).")

    # Pre-calculate all-pairs shortest paths for efficiency
    shortest_paths = dict(nx.all_pairs_shortest_path_length(G))

    print("\nStep 3: Calculating the Wiener Index (W)...")
    # The Wiener index is the sum of all shortest path distances between node pairs.
    # We sum over all ordered pairs (A->B and B->A) and then divide by 2.
    wiener_index = sum(sum(d.values()) for d in shortest_paths.values()) // 2
    
    print("\nStep 4: Calculating the Szeged Index (Sz)...")
    szeged_index = 0
    # Iterate over each edge (bond) in the graph
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        # For each edge, partition the vertices of the graph
        for w in G.nodes():
            dist_u_w = shortest_paths[u][w]
            dist_v_w = shortest_paths[v][w]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v

    print("\nStep 5: Calculating the final ratio...")
    if wiener_index == 0:
        ratio = 0
    else:
        ratio = szeged_index / wiener_index
        
    print("\n--- Calculation Complete ---")
    print(f"Molecule: Perylene-3-thiol (C20H11SH)")
    print(f"Wiener Index (W): {wiener_index}")
    print(f"Szeged Index (Sz): {szeged_index}")
    print("\nSzeged/Wiener Index Ratio:")
    print(f"{szeged_index} / {wiener_index} = {ratio}")

if __name__ == '__main__':
    # To run this code, you need to have rdkit and networkx installed.
    # You can install them using pip:
    # pip install rdkit-pypi networkx
    try:
        solve_szeged_wiener_ratio()
    except ImportError:
        print("Error: RDKit or NetworkX library not found.")
        print("Please install them using: pip install rdkit-pypi networkx")
