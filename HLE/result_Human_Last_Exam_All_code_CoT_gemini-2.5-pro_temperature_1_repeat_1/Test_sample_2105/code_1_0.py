import networkx as nx
from rdkit import Chem

def calculate_indices_and_ratio():
    """
    Calculates the Wiener and Szeged indices for perylene-3-thiol and their ratio.
    """
    # Step 1: Define the molecule and create the graph.
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    # We use its SMILES string to generate the structure.
    # SMILES for perylene-3-thiol: c1cc(S)c2c3ccc4cccc5ccc(c1c35)c24
    smiles_string = 'c1cc(S)c2c3ccc4cccc5ccc(c1c35)c24'

    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Could not create molecule from SMILES string.")
        return

    # Add hydrogen atoms to the graph as specified
    mol_with_hydrogens = Chem.AddHs(mol)

    # Create a NetworkX graph from the molecule
    G = nx.Graph()
    for atom in mol_with_hydrogens.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hydrogens.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    num_atoms = G.number_of_nodes()
    
    # Pre-calculate all-pairs shortest path lengths for efficiency
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(G.nodes())

    # Step 2: Calculate Wiener Index (W)
    # W = 0.5 * sum of all shortest path distances between any two vertices
    wiener_index = 0
    for start_node in nodes:
        for end_node in nodes:
            wiener_index += path_lengths[start_node][end_node]
    wiener_index //= 2  # Each path is counted twice, so divide by 2

    # Step 3: Calculate Szeged Index (Sz)
    # Sz = sum over all edges e=(u,v) of n_u(e) * n_v(e)
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in nodes:
            dist_u_w = path_lengths[u][w]
            dist_v_w = path_lengths[v][w]
            if dist_u_w < dist_v_w:
                n_u += 1
            elif dist_v_w < dist_u_w:
                n_v += 1
        szeged_index += n_u * n_v

    # Step 4: Compute and print the ratio
    if wiener_index > 0:
        ratio = szeged_index / wiener_index
        print(f"Analysis of Perylene-3-thiol (including all {num_atoms} atoms):")
        print("-" * 40)
        print(f"Calculated Wiener Index (W): {wiener_index}")
        print(f"Calculated Szeged Index (Sz): {szeged_index}")
        print("-" * 40)
        print("Final Ratio Calculation:")
        print(f"Sz / W = {szeged_index} / {wiener_index} = {ratio}")
    else:
        print("Cannot calculate ratio because Wiener index is zero.")

if __name__ == '__main__':
    calculate_indices_and_ratio()