import networkx as nx
import numpy as np

def analyze_graph(G, graph_name):
    """
    Analyzes a graph to see if it satisfies the given condition and checks the options.
    """
    print(f"--- Analyzing graph: {graph_name} ---")

    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    # In networkx, components are returned as node sets
    components = [G.subgraph(c) for c in nx.connected_components(G)]
    k = len(components)
    
    k_bip = 0
    for comp in components:
        if nx.is_bipartite(comp):
            k_bip += 1
            
    # Construct unoriented incidence matrix B
    # Note: networkx returns a sparse matrix
    B = nx.incidence_matrix(G, oriented=False).toarray()
    
    # Calculate nullity of B^T * B
    # nullity = num_columns - rank
    BTB = B.T @ B
    rank_BTB = np.linalg.matrix_rank(BTB)
    nullity_BTB = BTB.shape[1] - rank_BTB
    
    # Calculate nullity from our derived formula
    nullity_formula = m - n + k_bip
    
    # Calculate lambda_n from the Laplacian
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = np.max(eigenvalues)

    # Check the given condition
    lhs = nullity_BTB
    rhs = lambda_n / 2
    
    print(f"Graph Properties: n={n}, m={m}, k={k}, k_bip={k_bip}")
    print(f"Calculated null(B^T B) from matrix rank: {lhs}")
    print(f"Calculated null(B^T B) from formula (m-n+k_bip): {nullity_formula}")
    print(f"Equation to check: {lhs} == {rhs} ?")

    # Use a small tolerance for float comparison
    if not np.isclose(lhs, rhs):
        print("Condition does not hold for this graph.\n")
        return
        
    print("Condition holds for this graph.\n")

    print("Checking answer choices:")
    # Choice C: The graph is connected
    print(f"C. Is graph connected? (k=1): {k == 1}")

    # Choice D: The graph has exactly lambda_n(G)/2 connected components
    print(f"D. Does k == lambda_n/2? {k} == {rhs} ? {np.isclose(k, rhs)}")

    # Choice B: The graph has at least lambda_n(G)/2 connected components
    print(f"B. Does k >= lambda_n/2? {k} >= {rhs} ? {k >= rhs or np.isclose(k,rhs)}")
    print("-" * 30 + "\n")


# Example 1: K4 (Complete graph on 4 vertices)
G1 = nx.complete_graph(4)
analyze_graph(G1, "K4")

# Example 2: Two disjoint squares (C4 U C4)
G2 = nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
analyze_graph(G2, "C4 U C4")
