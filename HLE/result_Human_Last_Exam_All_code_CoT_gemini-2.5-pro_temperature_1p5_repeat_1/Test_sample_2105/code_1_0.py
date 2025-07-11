import networkx as nx

def solve():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    This script requires the rdkit and networkx libraries.
    You can install them using pip:
    pip install rdkit-pypi networkx
    """
    try:
        from rdkit import Chem
    except ImportError:
        print("RDKit library not found.")
        print("Please install it using: pip install rdkit-pypi")
        return

    # 1. Define the molecule, which is the major reduction product of
    #    di(perylene-3-yl) disulfide -> perylene-3-thiol.
    # We use its SMILES representation.
    smiles_perylene_thiol = "c1cc(S)c2c(c1)c3c4ccccc4c5cccc(c2)c35"

    # 2. Create a molecular object from SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles_perylene_thiol)
    mol_with_h = Chem.AddHs(mol)

    # 3. Convert the rdkit molecule to a networkx graph.
    #    Atoms are nodes, bonds are edges.
    G = nx.Graph()
    for atom in mol_with_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Pre-calculate all-pairs shortest paths for efficiency.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # 4. Calculate the Wiener Index (W).
    #    W = 1/2 * sum(d(i,j)) for all pairs of vertices (i, j).
    wiener_index = sum(val for dest in path_lengths.values() for val in dest.values()) / 2

    # 5. Calculate the Szeged Index (Sz).
    #    Sz = sum_{edges e=(u,v)} n_u(e) * n_v(e).
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

    # Ensure indices are integers for clear output
    wiener_index = int(wiener_index)
    szeged_index = int(szeged_index)

    # 6. Calculate the final ratio.
    if wiener_index == 0:
        ratio = 0
    else:
        ratio = szeged_index / wiener_index

    # Print the results in the required equation format.
    print(f"Molecule: Perylene-3-thiol (C20H12S)")
    print(f"Total atoms (nodes): {G.number_of_nodes()}")
    print(f"Total bonds (edges): {G.number_of_edges()}")
    print("-" * 30)
    print(f"Wiener Index (W): {wiener_index}")
    print(f"Szeged Index (Sz): {szeged_index}")
    print("-" * 30)
    print("Szeged/Wiener Ratio Calculation:")
    print(f"{szeged_index} / {wiener_index} = {ratio}")

if __name__ == "__main__":
    solve()