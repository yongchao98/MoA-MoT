import networkx as nx
from rdkit import Chem

# Memoization cache for the Hosoya Z index calculation
hosoya_memo = {}

def get_graph_from_smiles(smiles: str) -> nx.Graph:
    """Converts a SMILES string to a NetworkX graph."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")
    
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        # Treat all bonds as single edges for topological indices
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return graph

def calculate_zagreb1_index(graph: nx.Graph) -> int:
    """Calculates the first Zagreb index (M1)."""
    m1 = 0
    for i in graph.nodes():
        m1 += graph.degree(i) ** 2
    return m1

def calculate_hosoya_z_index(graph: nx.Graph) -> int:
    """
    Calculates the Hosoya Z index (Z) of a graph using a recursive formula
    with memoization. Z is the total number of matchings in the graph.
    """
    # The key for memoization is a canonical representation of the graph's edges.
    # A frozenset of sorted edge tuples is used.
    edge_set = frozenset(tuple(sorted(e)) for e in graph.edges())

    if not edge_set:
        return 1
    
    if edge_set in hosoya_memo:
        return hosoya_memo[edge_set]

    # Recursive step: Z(G) = Z(G-e) + Z(G-{u,v})
    # Pick an arbitrary edge e = (u,v)
    u, v = next(iter(edge_set))
    
    # Create G-e (graph with edge e removed)
    graph_minus_e = graph.copy()
    graph_minus_e.remove_edge(u, v)
    
    # Create G-{u,v} (graph with nodes u, v and incident edges removed)
    graph_minus_uv = graph.copy()
    graph_minus_uv.remove_nodes_from([u, v])

    # The result is the sum of the Hosoya indices of the two smaller graphs.
    result = calculate_hosoya_z_index(graph_minus_e) + calculate_hosoya_z_index(graph_minus_uv)
    hosoya_memo[edge_set] = result
    
    return result

def solve_task():
    """
    Solves the user's request by identifying the target molecule and calculating the required ratio.
    """
    # The target molecule identified is Deoxycytidine.
    # SMILES representation for Deoxycytidine (dC).
    smiles_dc = 'NC1=NC(=O)N=CN1C2CC(O)C(CO)O2'

    try:
        # Create the molecular graph
        molecular_graph = get_graph_from_smiles(smiles_dc)
        
        # Calculate Zagreb(1) index
        zagreb1_index = calculate_zagreb1_index(molecular_graph)
        
        # Clear memoization cache before new calculation
        global hosoya_memo
        hosoya_memo.clear()
        
        # Calculate Hosoya Z index
        hosoya_z_index = calculate_hosoya_z_index(molecular_graph)
        
        # Calculate the final ratio
        numerator_factor = 2
        final_ratio = (numerator_factor * hosoya_z_index) / zagreb1_index
        
        # Print the final equation as requested
        print(f"Target Molecule: Deoxycytidine (dC)")
        print(f"Zagreb(1) Index (M1): {zagreb1_index}")
        print(f"Hosoya Z Index (Z): {hosoya_z_index}")
        print("\nFinal Calculation:")
        print(f"({numerator_factor} * {hosoya_z_index}) / {zagreb1_index} = {final_ratio:.4f}")

    except (ValueError, ImportError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have rdkit and networkx installed (`pip install rdkit-pypi networkx`).")

if __name__ == "__main__":
    solve_task()
<<<28.3333>>>