import rdkit
from rdkit import Chem
import networkx as nx

def solve_molecular_indices():
    """
    Calculates the Wiener and Szeged indices for perylene-3-thiol and their ratio.
    """
    # Step 1 & 2: Define the molecule (perylene-3-thiol) and create its graph representation
    # The SMILES string for perylene-3-thiol is Sc1cc2c3c(ccc2cc1)c4cccc5c4c3ccc5
    smiles = "Sc1cc2c3c(ccc2cc1)c4cccc5c4c3ccc5"
    mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogen atoms to the molecular graph, as required
    mol = Chem.AddHs(mol)
    
    # Create a NetworkX graph from the RDKit molecule object
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Calculate all-pairs shortest paths
    # The result is a dictionary where dist[u][v] gives the shortest path length
    dist = dict(nx.all_pairs_shortest_path_length(G))

    # Step 4: Calculate the Wiener Index (W)
    # Sum of distances between all unique pairs of nodes
    wiener_index = 0
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            wiener_index += dist[u][v]

    # Step 5: Calculate the Szeged Index (Sz)
    # Sum of n_u * n_v for each edge (u, v)
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in G.nodes():
            dist_wu = dist[w][u]
            dist_wv = dist[w][v]
            if dist_wu < dist_wv:
                n_u += 1
            elif dist_wv < dist_wu:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 6: Compute and print the ratio
    if wiener_index == 0:
        ratio = 0
    else:
        ratio = szeged_index / wiener_index

    print(f"Molecule: Perylene-3-thiol (C20H12S)")
    print(f"Total atoms (including H): {mol.GetNumAtoms()}")
    print("-" * 30)
    print(f"Szeged Index (Sz): {szeged_index}")
    print(f"Wiener Index (W): {wiener_index}")
    print("-" * 30)
    print(f"Szeged Index / Wiener Index = {szeged_index} / {wiener_index} = {ratio}")

# Execute the function
solve_molecular_indices()