# First, ensure you have rdkit and networkx installed:
# pip install rdkit networkx

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem

def solve_indices_ratio():
    """
    This function calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    
    # Step 1: Define the molecule (perylene-3-thiol) and add hydrogens.
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    smiles = "Sc1cc2c3ccc4cccc5c4c(c3cc1)ccc5"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Could not create molecule from SMILES string.")
        return
    mol_h = Chem.AddHs(mol)

    # Step 2: Convert the RDKit molecule object to a NetworkX graph.
    graph = nx.Graph()
    for atom in mol_h.GetAtoms():
        graph.add_node(atom.GetIdx())
    for bond in mol_h.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Step 3: Calculate the Wiener Index (W).
    # The Wiener index is the sum of shortest path lengths between all pairs of nodes.
    try:
        wiener_index = nx.wiener_index(graph)
    except nx.NetworkXError as e:
        print(f"Error calculating Wiener Index: {e}")
        return

    # Step 4: Calculate the Szeged Index (Sz).
    def calculate_szeged_index(G):
        """Calculates the Szeged index of a graph G."""
        sz = 0
        
        # Pre-computing all-pairs shortest paths is more efficient.
        try:
            all_pairs_sp_length = dict(nx.all_pairs_shortest_path_length(G))
        except nx.NetworkXError as e:
            print(f"Error computing all-pairs shortest path: {e}")
            return None

        for u, v in G.edges():
            n_u = 0
            n_v = 0
            for w in G.nodes():
                dist_u = all_pairs_sp_length[w][u]
                dist_v = all_pairs_sp_length[w][v]
                
                if dist_u < dist_v:
                    n_u += 1
                elif dist_v < dist_u:
                    n_v += 1
            sz += n_u * n_v
        return sz

    szeged_index = calculate_szeged_index(graph)
    if szeged_index is None:
        return

    # Step 5: Calculate the ratio and print the final result.
    if wiener_index == 0:
        print("Wiener index is zero, ratio cannot be calculated.")
        return
        
    ratio = szeged_index / wiener_index
    
    print("Molecule: Perylene-3-thiol (including H atoms)")
    print(f"Number of atoms (vertices): {graph.number_of_nodes()}")
    print(f"Number of bonds (edges): {graph.number_of_edges()}")
    print("-" * 30)
    print("Szeged Index (Sz):")
    print(szeged_index)
    print("Wiener Index (W):")
    print(wiener_index)
    print("-" * 30)
    print("Szeged/Wiener Ratio (Sz / W):")
    # The final equation as requested
    print(f"{szeged_index} / {wiener_index} = {ratio}")
    
    # Store the final numerical answer for the platform.
    global final_answer
    final_answer = ratio

# Execute the function
final_answer = 0.0
solve_indices_ratio()
# The 'final_answer' variable holds the result in case it needs to be accessed programmatically.
# For example: print(f'Final numeric answer: {final_answer}')
