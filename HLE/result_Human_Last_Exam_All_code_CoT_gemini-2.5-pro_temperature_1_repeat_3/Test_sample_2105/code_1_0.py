import rdkit
from rdkit import Chem
import networkx as nx

def calculate_molecular_indices():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1: Define the molecule (perylene-3-thiol) using its SMILES string
    # and add explicit hydrogen atoms.
    smiles = "SC1=CC2=C3C=CC4=CC=CC5=C4C3=C(C=C2)C2=C(C=C1)C=CC5=C2"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Error: Could not create molecule from SMILES string.")
            return
        mol = Chem.AddHs(mol)
    except ImportError:
        print("Error: RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
        return


    # Step 2: Convert the RDKit molecule object into a NetworkX graph.
    G = nx.Graph()
    G.add_nodes_from(range(mol.GetNumAtoms()))
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Calculate the shortest path lengths between all pairs of atoms.
    # This is a prerequisite for both Wiener and Szeged indices.
    try:
        path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    except ImportError:
        print("Error: NetworkX is not installed. Please install it using 'pip install networkx'")
        return

    # Step 4: Calculate the Wiener Index (W).
    # The sum of distances between all unordered pairs of vertices.
    wiener_index = 0
    for source_node in path_lengths:
        for dest_node in path_lengths[source_node]:
            wiener_index += path_lengths[source_node][dest_node]
    wiener_index //= 2  # Divide by 2 for unordered pairs.

    # Step 5: Calculate the Szeged Index (Sz).
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

    # Step 6: Calculate the ratio and print the results in the required format.
    if wiener_index == 0:
        ratio = 0.0
    else:
        ratio = szeged_index / wiener_index

    print(f"Analysis for the major reduction product of di(perylene-3-yl) disulfide (Perylene-3-thiol):")
    print(f"Total number of atoms (vertices), including H: {mol.GetNumAtoms()}")
    print("-" * 30)
    print(f"Szeged Index (Sz): {szeged_index}")
    print(f"Wiener Index (W): {wiener_index}")
    print("-" * 30)
    print(f"Sz / W = {szeged_index} / {wiener_index} = {ratio}")

if __name__ == '__main__':
    calculate_molecular_indices()