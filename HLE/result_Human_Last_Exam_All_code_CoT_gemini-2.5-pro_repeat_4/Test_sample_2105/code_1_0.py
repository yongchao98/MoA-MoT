import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import networkx as nx

def calculate_indices():
    """
    Calculates the Wiener and Szeged indices for perylene-3-thiol.
    """
    # Step 1: Define and create the molecule perylene-3-thiol.
    # The SMILES string for perylene-3-thiol is c1ccc2c(c1)c3cc4c5c(cccc5)c(c4c(S)c3)c2
    smiles = "c1ccc2c(c1)c3cc4c5c(cccc5)c(c4c(S)c3)c2"
    mol = Chem.MolFromSmiles(smiles)
    
    # Step 2: Add explicit hydrogens as per the problem statement "including H".
    mol = Chem.AddHs(mol)
    
    # Create a NetworkX graph from the RDKit molecule object.
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Pre-calculate all-pairs shortest paths, which is needed for both indices.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))

    # Step 3: Calculate the Wiener Index (W).
    # W is the sum of shortest paths between all pairs of vertices, divided by 2.
    wiener_index = 0
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            wiener_index += path_lengths[nodes[i]][nodes[j]]

    # Step 4: Calculate the Szeged Index (Sz).
    # Sz = sum over all edges(u,v) of n_u * n_v
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
        
    # Step 5: Compute the ratio and print the results.
    ratio = szeged_index / wiener_index

    print(f"The major reduction product is perylene-3-thiol.")
    print(f"Molecular Formula: {rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"Total number of atoms (including H): {mol.GetNumAtoms()}")
    print("-" * 30)
    print(f"Calculated Wiener Index (W): {int(wiener_index)}")
    print(f"Calculated Szeged Index (Sz): {szeged_index}")
    print("-" * 30)
    print("Final Ratio Calculation:")
    print(f"{szeged_index} / {int(wiener_index)} = {ratio}")

if __name__ == '__main__':
    calculate_indices()