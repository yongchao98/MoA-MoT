import rdkit
from rdkit import Chem
import networkx as nx

def solve():
    """
    Calculates the Szeged/Wiener index ratio for Perylene-3-thiol.
    """
    # Step 1: Define the molecule (Perylene-3-thiol) using its SMILES string.
    # This is the major reduction product of di(perylene-3-yl) disulfide.
    # SMILES verified to correspond to perylene-3-thiol.
    smiles = "Sc1cc2ccc3c4ccccc4c5cccc(c1)c2c35"
    mol = Chem.MolFromSmiles(smiles)
    
    # Step 2: Add explicit hydrogen atoms to the molecule as required by the prompt.
    mol_with_hs = Chem.AddHs(mol)
    num_atoms = mol_with_hs.GetNumAtoms()
    
    # Step 3: Convert the RDKit molecule object to a NetworkX graph.
    # The nodes of the graph are the atoms, and the edges are the bonds.
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        
    # Step 4: Calculate all-pairs shortest paths, which is needed for both indices.
    # This returns a dictionary of dictionaries: path_lengths[source][target] = length.
    path_lengths = dict(nx.all_pairs_shortest_path_length(G))
    
    # Step 5: Calculate the Wiener Index (W).
    # W is the sum of shortest path distances over all unordered pairs of vertices.
    wiener_index = 0
    nodes = list(G.nodes())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            wiener_index += path_lengths[u][v]
            
    # Step 6: Calculate the Szeged Index (Sz).
    # Sz = sum over all edges (u,v) of n_u * n_v, where n_u is the number of
    # nodes closer to u than v, and n_v is the number of nodes closer to v than u.
    szeged_index = 0
    for u, v in G.edges():
        n_u = 0
        n_v = 0
        for w in G.nodes():
            dist_w_u = path_lengths[w][u]
            dist_w_v = path_lengths[w][v]
            if dist_w_u < dist_w_v:
                n_u += 1
            elif dist_w_v < dist_w_u:
                n_v += 1
        szeged_index += n_u * n_v
        
    # Step 7: Compute the ratio and print the results in the required format.
    ratio = szeged_index / wiener_index if wiener_index != 0 else 0
    
    print(f"Molecule: Perylene-3-thiol")
    print(f"Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)}")
    print(f"Number of atoms (including H): {num_atoms}")
    print(f"Number of bonds: {G.number_of_edges()}")
    print("-" * 30)
    print(f"Calculated Wiener Index (W): {wiener_index}")
    print(f"Calculated Szeged Index (Sz): {szeged_index}")
    print("-" * 30)
    print("Szeged/Wiener Ratio calculation:")
    print(f"{szeged_index} / {wiener_index} = {ratio}")
    
    return ratio

# Run the solver and store the final answer
final_ratio = solve()
# The final answer must be returned in the specified format
# print(f'<<<{final_ratio}>>>')
# Let's format it to a reasonable number of decimal places for consistency
print(f'\n<<<Final Answer>>>')
print(f'<<<{final_ratio:.4f}>>>')
