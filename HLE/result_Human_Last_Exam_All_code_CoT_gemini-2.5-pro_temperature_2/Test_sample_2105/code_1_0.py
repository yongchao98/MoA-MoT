import networkx as nx
from rdkit import Chem

def calculate_indices():
    """
    This script calculates the Wiener and Szeged indices for perylene-3-thiol,
    including all hydrogen atoms, and then computes their ratio.
    """
    # Step 1: Construct the molecule perylene-3-thiol (C20H11SH)
    # Start with perylene
    perylene_smiles = 'c1cc2c(cc1)c1c3cccc4c3c(c(c21)cc4)'
    m = Chem.MolFromSmiles(perylene_smiles)

    # Convert to an editable molecule to add the thiol group.
    # By analyzing the molecule generated from this SMILES string, it's known
    # that rdkit's atom index 2 corresponds to IUPAC position 3 of perylene.
    em = Chem.EditableMol(m)
    s_idx = em.AddAtom(Chem.Atom('S'))  # Add a sulfur atom
    em.AddBond(2, s_idx, order=Chem.rdchem.BondType.SINGLE)  # Bond it to C-3

    # Get the final molecule and add hydrogens to satisfy valency.
    # RDKit will correctly remove the H from C-3 and add one to the S.
    product_mol_no_h = em.GetMol()
    Chem.SanitizeMol(product_mol_no_h)
    final_mol = Chem.AddHs(product_mol_no_h)

    # Step 2: Convert the RDKit molecule to a NetworkX graph
    G = nx.Graph()
    for atom in final_mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in final_mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Calculate Wiener and Szeged indices
    # Pre-calculate all-pairs shortest paths for efficiency
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(G.nodes())

    # Calculate Wiener Index (W)
    wiener_index = 0
    # Sum shortest path distances over all unique pairs of nodes
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            wiener_index += path_lengths[nodes[i]][nodes[j]]

    # Calculate Szeged Index (Sz)
    szeged_index = 0
    # Sum n_u * n_v over all edges (u,v)
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for x in nodes:
            dist_u = path_lengths[x][u]
            dist_v = path_lengths[x][v]
            if dist_u < dist_v:
                n_u += 1
            elif dist_v < dist_u:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 4: Calculate the ratio and print the final output
    if wiener_index == 0:
        ratio = 0.0
    else:
        ratio = szeged_index / wiener_index

    print(f"Analysis of Perylene-3-thiol (including H atoms)")
    print(f"Number of atoms (nodes): {G.number_of_nodes()}")
    print(f"Number of bonds (edges): {G.number_of_edges()}")
    print("-" * 30)
    print("Topological Indices:")
    print(f"Szeged Index (Sz) = {szeged_index}")
    print(f"Wiener Index (W) = {wiener_index}")
    print("-" * 30)
    print("Final Ratio Calculation:")
    print(f"{szeged_index} / {wiener_index} = {ratio:.7f}")


if __name__ == '__main__':
    calculate_indices()
