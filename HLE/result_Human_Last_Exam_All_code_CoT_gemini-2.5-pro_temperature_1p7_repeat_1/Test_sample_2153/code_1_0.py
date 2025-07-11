# You may need to install the required libraries:
# pip install rdkit-pypi networkx

from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors, rdMolDescriptors
import networkx as nx
import sys

# Increase recursion limit for the Hosoya Z calculation, as it can be intensive.
sys.setrecursionlimit(2000)

# Memoization cache for the recursive Hosoya Z function to improve performance.
hosoya_memo = {}

def calculate_hosoya_z(graph):
    """
    Calculates the Hosoya Z index of a graph using a recursive algorithm with memoization.
    The Hosoya Z index is the total number of matchings (sets of non-adjacent edges).
    Z(G) = Z(G-e) + Z(G-V(e))
    """
    # Create a canonical, hashable representation of the graph's edges for memoization.
    graph_edges_tuple = tuple(sorted(map(tuple, sorted(graph.edges()))))
    
    if graph_edges_tuple in hosoya_memo:
        return hosoya_memo[graph_edges_tuple]

    if graph.number_of_edges() == 0:
        return 1

    # Select an edge to start the recursion
    edge = list(graph.edges())[0]
    u, v = edge

    # Term 1: Graph without the selected edge 'e'
    graph_minus_edge = graph.copy()
    graph_minus_edge.remove_edge(u, v)
    term1 = calculate_hosoya_z(graph_minus_edge)

    # Term 2: Graph without the vertices of the selected edge 'e'
    graph_minus_vertices = graph.copy()
    graph_minus_vertices.remove_nodes_from([u, v])
    term2 = calculate_hosoya_z(graph_minus_vertices)
    
    result = term1 + term2
    hosoya_memo[graph_edges_tuple] = result
    return result

def solve_chemistry_problem():
    """
    Main function to orchestrate the problem-solving steps.
    """
    # Step 1: Identify the BCKDH substrate with median Bertz complexity.
    bcaas = {
        'Leucine': 'CC(C)CC(C(=O)O)N',
        'Isoleucine': 'CCC(C)C(C(=O)O)N',
        'Valine': 'CC(C)C(C(=O)O)N'
    }
    bcaas_bertz = {}
    for name, smiles in bcaas.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            bcaas_bertz[name] = Descriptors.BertzCT(mol)
    
    # Sort by Bertz complexity to find the median
    sorted_bcaas = sorted(bcaas_bertz.items(), key=lambda item: item[1])
    median_bcaa_name = sorted_bcaas[1][0]
    
    # Step 2: Identify the nucleoside with the closest Balaban J index.
    median_bcaa_mol = Chem.MolFromSmiles(bcaas[median_bcaa_name])
    target_balaban_j = GraphDescriptors.BalabanJ(median_bcaa_mol)

    nucleosides = {
        'Adenosine': 'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)CO)O)O)N',
        'Guanosine': 'c1nc2c(n1C3C(C(C(O3)CO)O)O)N=C(NC2=O)N',
        'Cytidine': 'c1cn(c(=O)nc1N)C2C(C(C(O2)CO)O)O',
        'Uridine': 'c1cn(c(=O)nc1=O)C2C(C(C(O2)CO)O)O'
    }
    
    balaban_differences = {}
    for name, smiles in nucleosides.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            j_index = GraphDescriptors.BalabanJ(mol)
            balaban_differences[name] = abs(j_index - target_balaban_j)

    target_substance_name = min(balaban_differences, key=balaban_differences.get)
    print(f"Target substance identified: {target_substance_name}")

    # Step 3 & 4: Calculate indices for the target substance and find the final ratio.
    target_smiles = nucleosides[target_substance_name]
    mol_heavy_atom = Chem.MolFromSmiles(target_smiles)

    # Calculate Zagreb(1) index on the heavy-atom graph
    zagreb1_index = rdMolDescriptors.CalcZagreb1(mol_heavy_atom)
    print(f"Zagreb(1) index: {int(zagreb1_index)}")

    # Calculate Hosoya Z index on the H-included graph
    mol_h_included = Chem.AddHs(mol_heavy_atom)
    
    # Convert RDKit molecule to NetworkX graph for Hosoya Z calculation
    nx_graph = nx.Graph()
    for atom in mol_h_included.GetAtoms():
        nx_graph.add_node(atom.GetIdx())
    for bond in mol_h_included.GetBonds():
        nx_graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    hosoya_z_index = calculate_hosoya_z(nx_graph)
    print(f"Hosoya Z index (H-included): {hosoya_z_index}")

    # Final calculation
    final_result = (2 * hosoya_z_index) / zagreb1_index
    
    print(f"Final calculation: (2 * {hosoya_z_index}) / {int(zagreb1_index)} = {final_result}")

if __name__ == "__main__":
    solve_chemistry_problem()
<<<84.5625>>>