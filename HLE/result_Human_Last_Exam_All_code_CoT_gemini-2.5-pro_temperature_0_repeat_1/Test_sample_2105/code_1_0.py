import networkx as nx
from rdkit import Chem

def solve():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1 & 2: Identify the molecule and construct its graph representation.
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    # SMILES for perylene-3-thiol.
    smi = 'C1=CC=C2C3=C4C=CC(=CC4=C(C=C2)C2=C1C=CC=C2)S'
    
    # Create molecule from SMILES and add hydrogen atoms.
    mol = Chem.MolFromSmiles(smi)
    mol_with_hs = Chem.AddHs(mol)
    
    # Create a networkx graph from the molecule.
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Pre-calculate all-pairs shortest paths for efficiency.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # Step 3: Calculate the Wiener Index (W).
    # The Wiener index is the sum of all shortest paths, divided by 2.
    wiener_index = 0
    for source_node in path_lengths:
        for target_node in path_lengths[source_node]:
            wiener_index += path_lengths[source_node][target_node]
    wiener_index //= 2

    # Step 4: Calculate the Szeged Index (Sz).
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

    # Step 5: Compute the ratio and print the results.
    ratio = szeged_index / wiener_index
    
    print("The molecule is perylene-3-thiol, the reduction product of di(perylene-3-yl) disulfide.")
    print(f"The molecular graph includes all {mol_with_hs.GetNumAtoms()} atoms (C, S, and H).")
    print("-" * 30)
    print(f"Calculated Wiener Index (W): {wiener_index}")
    print(f"Calculated Szeged Index (Sz): {szeged_index}")
    print("-" * 30)
    print(f"Szeged/Wiener Ratio = {szeged_index} / {wiener_index} = {ratio}")

solve()
<<<3.456937799043062>>>