import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

# It might be necessary to increase recursion limit for Hosoya index calculation
sys.setrecursionlimit(2000)

# Dictionary for caching Hosoya Z index results to avoid re-computation
memo = {}

def hosoya_z_index(graph):
    """
    Calculates the Hosoya Z index (number of matchings) of a graph using a recursive algorithm.
    Memoization is used to cache results for subgraphs to ensure efficiency.
    """
    # Use a frozenset of edges as a unique, hashable identifier for the graph
    graph_id = frozenset(graph.edges())
    if graph_id in memo:
        return memo[graph_id]

    # Base case: A graph with no edges has one matching (the empty set)
    if graph.number_of_edges() == 0:
        return 1

    # Pick an edge to start the recursion
    u, v = next(iter(graph.edges()))

    # Z(G) = Z(G-e) + Z(G-N(e))
    # Term 1: Calculate Z for the graph without the chosen edge 'e'
    graph_minus_edge = graph.copy()
    graph_minus_edge.remove_edge(u, v)
    term1 = hosoya_z_index(graph_minus_edge)

    # Term 2: Calculate Z for the graph without the vertices of 'e' and their incident edges
    graph_minus_nodes = graph.copy()
    graph_minus_nodes.remove_nodes_from([u, v])
    term2 = hosoya_z_index(graph_minus_nodes)

    result = term1 + term2
    memo[graph_id] = result
    return result

def solve_chemoinformatics_puzzle():
    """
    Main function to execute all steps of the puzzle.
    """
    # Step 1: Define SMILES for all relevant molecules
    bcaa_smiles = {
        'Valine': 'CC(C)C(N)C(=O)O',
        'Leucine': 'CC(C)CC(N)C(=O)O',
        'Isoleucine': 'CCC(C)C(N)C(=O)O'
    }

    nucleoside_smiles = {
        'Cytidine': 'C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O',
        'Uridine': 'C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O',
        'Deoxyadenosine': 'NC1=NC=NC2=C1N=CN2C3CC(O)C(CO)O3',
        'Deoxyguanosine': 'NC1=NC2=C(N=C1O)N=CN2C3CC(O)C(CO)O3'
    }

    # Step 2: Determine the BCKDH-related substrate with median Bertz complexity
    bcaa_data = []
    for name, smiles in bcaa_smiles.items():
        mol = Chem.MolFromSmiles(smiles)
        bertz_ct = GraphDescriptors.BertzCT(mol)
        bcaa_data.append({'name': name, 'bertz': bertz_ct, 'mol': mol})

    # Sort by Bertz complexity to find the median
    bcaa_data.sort(key=lambda x: x['bertz'])
    median_molecule_data = bcaa_data[1] # The middle element in a list of 3
    target_name = median_molecule_data['name']
    
    # Step 3: Calculate the target Balaban J index for the median molecule
    target_j_index = GraphDescriptors.BalabanJ(median_molecule_data['mol'])
    print(f"The molecule with median Bertz complexity is {target_name}.")
    print(f"Its Balaban J index is: {target_j_index:.4f}\n")

    # Step 4: Identify the nucleoside with the closest Balaban J index
    best_match_name = ''
    min_diff = float('inf')
    final_mol_smiles = ''

    print("Comparing Balaban J indices of nucleosides to the target...")
    for name, smiles in nucleoside_smiles.items():
        mol = Chem.MolFromSmiles(smiles)
        j_index = GraphDescriptors.BalabanJ(mol)
        diff = abs(j_index - target_j_index)
        print(f"- {name}: J = {j_index:.4f}, Difference = {diff:.4f}")
        if diff < min_diff:
            min_diff = diff
            best_match_name = name
            final_mol_smiles = smiles
    
    print(f"\nThe substance with the most similar Balaban J index is: {best_match_name}\n")

    # Step 5: Perform graph-theoretic calculations for the final molecule
    final_mol = Chem.MolFromSmiles(final_mol_smiles)

    # Convert to a NetworkX graph (hydrogen-depleted)
    G = nx.Graph()
    for atom in final_mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in final_mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Calculate Zagreb(1) index
    degrees = [val for _, val in G.degree()]
    zagreb_m1 = sum(d**2 for d in degrees)
    
    # Calculate Hosoya Z index
    hosoya_z = hosoya_z_index(G)
    
    # Step 6: Calculate and display the final ratio
    final_ratio = (2 * hosoya_z) / zagreb_m1
    
    print("Final Calculation:")
    print(f"Twice Hosoya Z to Zagreb(1) Index Ratio for {best_match_name}")
    print(f"Hosoya Z Index (Z) = {hosoya_z}")
    print(f"Zagreb(1) Index (M1) = {zagreb_m1}")
    print(f"Ratio = (2 * Z) / M1 = (2 * {hosoya_z}) / {zagreb_m1} = {final_ratio:.4f}")

solve_chemoinformatics_puzzle()
<<<27.2683>>>