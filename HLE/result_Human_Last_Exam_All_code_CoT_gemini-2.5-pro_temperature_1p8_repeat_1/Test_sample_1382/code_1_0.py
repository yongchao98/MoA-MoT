import numpy as np
import networkx as nx

def analyze_graph(G, name):
    """
    Analyzes a given graph to check the property null(B^T B) = lambda_n / 2
    and evaluates the given answer choices.
    """
    print(f"--- Analyzing graph: {name} ---")
    
    # Graph properties
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    # The number of connected components is c.
    # The number of bipartite components is c_b.
    is_bipartite_list = [nx.is_bipartite(G.subgraph(c)) for c in nx.connected_components(G)]
    c_b = sum(is_bipartite_list)

    # Calculate null(B^T B) = m - n + c_b
    nullity_BtB = m - n + c_b
    print(f"Graph properties: n={n}, m={m}, bipartite components c_b={c_b}")
    print(f"Calculated null(B^T B) = m - n + c_b = {m} - {n} + {c_b} = {nullity_BtB}")

    # Calculate lambda_n(G) / 2
    # Check if graph has edges, otherwise Laplacian is zero matrix
    if m > 0:
        L = nx.laplacian_matrix(G).toarray()
        eigenvalues = np.linalg.eigvalsh(L)
        lambda_n = np.max(eigenvalues)
    else: # No edges
        lambda_n = 0
    
    lambda_n_div_2 = lambda_n / 2
    print(f"Largest Laplacian eigenvalue lambda_n = {lambda_n:.4f}")
    print(f"lambda_n / 2 = {lambda_n_div_2:.4f}")

    # Check the condition from the problem
    if np.isclose(nullity_BtB, lambda_n_div_2):
        print("\nThe condition null(B^T B) = lambda_n / 2 HOLDS for this graph.")
    else:
        print("\nThe condition null(B^T B) = lambda_n / 2 DOES NOT HOLD for this graph.")
        return

    # Now evaluate the answer choices for this graph
    print("\n--- Evaluating Answer Choices ---")
    
    c = nx.number_connected_components(G)

    # A. If you drop lambda_n/2 edges, there will be at least two nodes with degree <= 1
    # This check is complex to generalize, but for K4, we drop 2 edges.
    # Removing two disjoint edges from K4 results in a C4, where all degrees are 2.
    # So A is false for K4. We will state this conclusion.
    print(f"A. For K4, lambda_n/2 = 2. Dropping 2 disjoint edges gives C4 (all degrees=2).")
    print("   So Choice A is False for K4.")

    # B. The graph has at least lambda_n(G)/2 connected components
    choice_b_holds = c >= lambda_n_div_2
    print(f"B. Is c >= lambda_n/2? -> Is {c} >= {lambda_n_div_2:.2f}? {choice_b_holds}")

    # C. The graph is connected
    choice_c_holds = (c == 1)
    print(f"C. Is the graph connected? c = {c}. {choice_c_holds}")
    
    # D. The graph has exactly lambda_n(G)/2 connected components
    choice_d_holds = np.isclose(c, lambda_n_div_2)
    print(f"D. Is c = lambda_n/2? -> Is {c} = {lambda_n_div_2:.2f}? {choice_d_holds}")


# Create K4, which we showed satisfies the condition
G1 = nx.complete_graph(4)
analyze_graph(G1, "K4 (Complete Graph on 4 nodes)")

# Based on the analysis, for K4, options A, B, and D are false. Option C is true.
# To show option C is not a general implication, we'd need the K4 U K2 example,
# which is a disconnected graph that also fulfills the condition.
# This proves that the property is not 'being connected'.

<<<E>>>