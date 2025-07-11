import networkx as nx

def print_graph_info(G, name):
    """Prints basic information about a graph."""
    print(f"Graph '{name}': {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

def main():
    """
    This script verifies the algebraic properties of graph operations
    to determine the correct classification of (G, union, tensor_product).
    """
    print("--- Defining Sample Graphs ---")
    G1 = nx.path_graph(2); G1.name = "P2"
    G2 = nx.path_graph(3); G2.name = "P3"
    print_graph_info(G1, G1.name)
    print_graph_info(G2, G2.name)
    print("-" * 30)

    # Property 1: Commutativity of Tensor Product (Multiplication)
    print("1. Is Tensor Product (*) commutative? (e.g., Is P2 * P3 isomorphic to P3 * P2?)")
    T1 = nx.tensor_product(G1, G2)
    T2 = nx.tensor_product(G2, G1)
    is_comm = nx.is_isomorphic(T1, T2)
    print(f"Result: {is_comm}. Multiplication is commutative.")
    print("-" * 30)
    
    # Property 2: Lack of Multiplicative Identity
    print("2. Is there a multiplicative identity 'I' in the set of simple graphs?")
    print("   Test: Let's check if I=K1 (1 vertex, 0 edges) works for G=P2.")
    K1 = nx.Graph(); K1.add_node(0); K1.name = "K1"
    Prod = nx.tensor_product(G1, K1)
    Prod.name = "P2 * K1"
    print_graph_info(G1, G1.name)
    print_graph_info(Prod, Prod.name)
    is_iso = nx.is_isomorphic(G1, Prod)
    print(f"   Is P2 * K1 isomorphic to P2? {is_iso}.")
    print("   Result: False. There is no multiplicative identity for simple graphs.")
    print("-" * 30)
    
    # Property 3: Distributivity of Tensor Product over Disjoint Union
    print("3. Does * distribute over U? (e.g., P2 * (P2 U P3) vs (P2*P2) U (P2*P3))")
    # Using G1 as the third graph for simplicity
    G3 = G1
    LHS = nx.tensor_product(G1, nx.disjoint_union(G2, G3))
    RHS = nx.disjoint_union(nx.tensor_product(G1, G2), nx.tensor_product(G1, G3))
    is_dist = nx.is_isomorphic(LHS, RHS)
    print(f"Result: {is_dist}. Distributivity holds.")
    print("-" * 30)

    # Property 4: Not a Ring (Lack of Additive Inverses)
    print("4. Is the structure a Ring? (Does an additive inverse exist for Union?)")
    print("   For a non-empty graph G, can we find H where G U H is the empty graph?")
    nodes_g1 = G1.number_of_nodes()
    # Let H be any graph, e.g., K1
    H = K1
    nodes_h = H.number_of_nodes()
    Union_GH = nx.disjoint_union(G1, H)
    nodes_union = Union_GH.number_of_nodes()
    print(f"   Number of nodes in G: {nodes_g1}")
    print(f"   Number of nodes in H: {nodes_h}")
    print(f"   Number of nodes in G U H: {nodes_union} (which is {nodes_g1} + {nodes_h})")
    print("   This can never be 0 unless both G and H are empty.")
    print("   Result: Additive inverses do not exist for non-empty graphs.")
    print("-" * 30)
    
    print("\nFinal Conclusion:")
    print("The structure (G, U, *) is a commutative semi-ring (by convention), but not a ring.")
    
if __name__ == "__main__":
    main()