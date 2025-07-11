import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

# It is known that the recursive calculation for Hosoya Z can reach Python's
# recursion depth limit for complex molecules. We increase it as a precaution.
sys.setrecursionlimit(2000)

# --- Part 1: Hosoya Z Index Calculation Function ---
# A memoization dictionary to store results for previously seen graphs
# to avoid re-computation, which is critical for performance.
_hosoya_memo = {}

def get_graph_key(graph):
    """Generates a canonical, hashable key for a graph based on its edges."""
    # A frozenset of frozensets of edge nodes is a canonical representation
    # for an undirected graph, suitable for dictionary keys.
    return frozenset(frozenset(edge) for edge in graph.edges())

def calculate_hosoya_z(graph):
    """
    Calculates the Hosoya Z index of a graph using a recursive formula with memoization.
    The Hosoya Z index is the total number of matchings (sets of non-adjacent edges) in a graph.
    The recurrence is Z(G) = Z(G-e) + Z(G-V(e)), where e is an edge and V(e) are its endpoints.
    """
    graph_key = get_graph_key(graph)
    if graph_key in _hosoya_memo:
        return _hosoya_memo[graph_key]

    if graph.number_of_edges() == 0:
        return 1
    
    # Pick an edge to recurse on. Any edge will do.
    edge = list(graph.edges())[0]
    u, v = edge

    # Create the graph G-e (remove the edge e)
    graph_minus_e = graph.copy()
    graph_minus_e.remove_edge(u, v)
    term1 = calculate_hosoya_z(graph_minus_e)

    # Create the graph G-V(e) (remove the vertices u and v)
    graph_minus_uv = graph.copy()
    graph_minus_uv.remove_nodes_from([u, v])
    term2 = calculate_hosoya_z(graph_minus_uv)

    result = term1 + term2
    _hosoya_memo[graph_key] = result
    return result

# --- Part 2: Main Logic to Solve the Problem ---

def solve_cheminformatics_problem():
    """
    This function follows the step-by-step plan to solve the user's request.
    """
    # Step 1: Identify the BCKDH complex substrate with median Bertz's complexity.
    # These are the keto-acids derived from branched-chain amino acids.
    bckdh_substrates = {
        'KIC': 'CC(C)CC(=O)C(=O)O',       # from Leucine
        'KMV': 'CCC(C)C(=O)C(=O)O',        # from Isoleucine
        'KIV': 'CC(C)C(=O)C(=O)O'         # from Valine
    }

    complexities = []
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        # BertzCT measures molecular complexity.
        bertz_ct = GraphDescriptors.BertzCT(mol)
        complexities.append({'name': name, 'smiles': smiles, 'bertz_ct': bertz_ct})
    
    # Sort by complexity to find the median.
    sorted_substrates = sorted(complexities, key=lambda x: x['bertz_ct'])
    median_substrate = sorted_substrates[1] # The middle one in a list of 3

    # Step 2: Calculate the Balaban J index for this median substrate to use as a reference.
    # Balaban J is a topological index based on distance sums.
    median_mol = Chem.MolFromSmiles(median_substrate['smiles'])
    target_balaban_j = GraphDescriptors.BalabanJ(median_mol)

    # Step 3: Identify the substance from Carrell's synthesis with a similar Balaban J index.
    # The substances are Inosine, Guanosine, and their N3-arabino isomers.
    carrell_substances = {
        'Inosine': 'O=c1[nH]c2ncn([C@H]3O[C@H](CO)[C@@H](O)[C@@H]3O)c2n1',
        'Guanosine': 'Nc1nc(=O)c2[nH]cnc2n1[C@H]3O[C@H](CO)[C@@H](O)[C@@H]3O',
        'iso-Inosine': 'O=c1[nH]cnc2c1ncn2[C@H]3O[C@H](CO)[C@H](O)[C@H]3O',
        'iso-Guanosine': 'Nc1nc(=O)c2[nH]cn(c2n1)[C@H]3O[C@H](CO)[C@H](O)[C@H]3O'
    }

    best_match = None
    min_diff = float('inf')

    for name, smiles in carrell_substances.items():
        mol = Chem.MolFromSmiles(smiles)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - target_balaban_j)
        if diff < min_diff:
            min_diff = diff
            best_match = {'name': name, 'smiles': smiles}

    # Step 4: Calculate the required indices for the identified substance (Guanosine).
    target_smiles = best_match['smiles']
    mol_target = Chem.MolFromSmiles(target_smiles)

    # Calculate Zagreb(1) index on the hydrogen-suppressed graph.
    zagreb1_index = GraphDescriptors.Zagreb1Index(mol_target)

    # Calculate Hosoya Z index on the hydrogen-included graph as requested.
    mol_target_h = Chem.AddHs(mol_target)
    
    # Convert the RDKit molecule with hydrogens to a NetworkX graph.
    G = nx.Graph()
    for atom in mol_target_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_target_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
    hosoya_z_index = calculate_hosoya_z(G)

    # Step 5: Calculate the final ratio and print the result.
    # Ratio = (2 * Hosoya Z) / Zagreb(1)
    final_ratio = (2 * hosoya_z_index) / zagreb1_index

    # Print the final result in the requested format.
    print(f"The target substance was identified as {best_match['name']}.")
    print(f"Its Zagreb(1) index is {zagreb1_index}.")
    print(f"Its H-included Hosoya Z index is {hosoya_z_index}.")
    print("\nThe final calculation is:")
    print(f"(2 * {hosoya_z_index}) / {zagreb1_index} = {final_ratio}")

# Execute the main function to find the solution.
solve_cheminformatics_problem()