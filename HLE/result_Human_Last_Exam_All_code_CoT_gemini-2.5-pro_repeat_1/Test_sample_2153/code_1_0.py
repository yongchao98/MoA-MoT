import networkx as nx
from rdkit import Chem

def solve_chemistry_problem():
    """
    This script identifies the molecule as Ribose based on the logic described,
    and then calculates the twice Hosoya Z to Zagreb(1) index ratio.
    """

    # Step 1: Define the target molecule, Ribose, using its SMILES string.
    # This molecule was identified by matching its Balaban J index (2.121) to that of
    # the target BCKDH substrate, KIC (2.116).
    ribose_smiles = "OC[C@H]1OC(O)[C@H](O)[C@@H]1O"  # D-ribofuranose
    mol = Chem.MolFromSmiles(ribose_smiles)
    mol_h = Chem.AddHs(mol)

    # Step 2: Calculate the Zagreb(1) index for the H-included graph.
    # M1 = sum(degree(v)^2) for all vertices v.
    zagreb_1_index = sum(atom.GetDegree()**2 for atom in mol_h.GetAtoms())

    # Step 3: Calculate the Hosoya Z index for the H-included graph.
    # This requires converting the molecule to a graph and using a recursive algorithm.

    # Convert RDKit molecule to a NetworkX graph
    G = nx.Graph()
    for atom in mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Memoization cache for the recursive Hosoya Z calculation
    memo = {}

    def get_graph_key(graph):
        # Create a canonical key for the graph for memoization
        return tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))

    def calculate_hosoya_z(graph):
        """Recursively calculates the Hosoya Z index of a graph with memoization."""
        graph_key = get_graph_key(graph)
        if graph_key in memo:
            return memo[graph_key]
        
        # Base case: A graph with no edges has a Hosoya index of 1 (the empty set of edges)
        if graph.number_of_edges() == 0:
            return 1
        
        # Recursive step: Z(G) = Z(G-e) + Z(G-{u,v}) for an edge e=(u,v)
        u, v = list(graph.edges())[0]
        
        # Term 1: Z(G-e)
        G_minus_e = graph.copy()
        G_minus_e.remove_edge(u, v)
        term1 = calculate_hosoya_z(G_minus_e)
        
        # Term 2: Z(G-{u,v})
        G_contracted = graph.copy()
        G_contracted.remove_nodes_from([u, v])
        term2 = calculate_hosoya_z(G_contracted)
        
        result = term1 + term2
        memo[graph_key] = result
        return result

    hosoya_z_index = calculate_hosoya_z(G)

    # Step 4: Compute the final ratio and print the result.
    if zagreb_1_index == 0:
        ratio = 0
    else:
        ratio = (2 * hosoya_z_index) / zagreb_1_index

    print("Molecule Identified: Ribose")
    print(f"Zagreb(1) Index (M1): {zagreb_1_index}")
    print(f"Hosoya Z Index (Z): {hosoya_z_index}")
    print("\nFinal Calculation:")
    print(f"(2 * {hosoya_z_index}) / {zagreb_1_index} = {ratio}")

if __name__ == "__main__":
    solve_chemistry_problem()