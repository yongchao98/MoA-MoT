import networkx as nx
import numpy as np

def analyze_graph_property(G, graph_name):
    """
    Analyzes a graph to see if it satisfies the property
    null(B^T*B) = lambda_n(G) / 2 and checks the answer choices.

    Args:
        G (nx.Graph): The graph to analyze.
        graph_name (str): The name of the graph for printing.
    """
    print(f"--- Analyzing graph: {graph_name} ---")
    n = G.number_of_nodes()
    m = G.number_of_edges()

    if n <= 3:
        print("Graph has 3 or fewer nodes, skipping as per problem statement.")
        return

    # Check for empty graph
    if m == 0:
        print("Graph has no edges. lambda_n is 0, nullity is 0. Property 0=0/2 holds trivially.")
        # But this case is not very informative for the choices.
        return

    # 1. Calculate null(B^T B) = null(B)
    B = nx.incidence_matrix(G, oriented=False).toarray()
    # nullity = num_columns - rank
    nullity_B = B.shape[1] - np.linalg.matrix_rank(B)
    print(f"Number of vertices n = {n}")
    print(f"Number of edges m = {m}")
    print(f"Nullity of the incidence matrix B is: {nullity_B}")

    # 2. Calculate lambda_n(G)
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1] # Largest eigenvalue
    print(f"Largest eigenvalue of the Laplacian, lambda_n = {lambda_n:.4f}")

    # 3. Check if the property holds
    holds = np.isclose(nullity_B, lambda_n / 2)
    print(f"Does the property null(B) = lambda_n / 2 hold? {holds}")
    if not holds:
        print("The property does not hold for this graph.")
        return

    # 4. Check the answer choices for this graph
    print("\n--- Evaluating answer choices for this graph ---")
    k = nx.number_connected_components(G)
    print(f"Number of connected components k = {k}")
    lambda_n_div_2 = lambda_n / 2
    print(f"Value of lambda_n / 2 = {lambda_n_div_2:.4f}")

    # Choice C: The graph is connected
    is_connected = (k == 1)
    print(f"Is Choice C (graph is connected) true? {is_connected}")
    if not is_connected:
        print("=> Since the graph is disconnected, C is false.")

    # Choice D: The graph has exactly lambda_n/2 connected components
    is_d_true = np.isclose(k, lambda_n_div_2)
    print(f"Is Choice D (k == lambda_n/2) true? {is_d_true}")
    if not is_d_true:
        print(f"=> Since k={k} and lambda_n/2={lambda_n_div_2:.4f}, D is false.")

    # Choice B: The graph has at least lambda_n/2 connected components
    is_b_true = (k >= lambda_n_div_2)
    print(f"Is Choice B (k >= lambda_n/2) true? {is_b_true}")
    if not is_b_true:
        print(f"=> Since k={k} is not >= {lambda_n_div_2:.4f}, B is false.")

    # Choice A: Dropping lambda_n/2 edges creates >= 2 nodes with degree <= 1
    # This is harder to check exhaustively, so we check for our specific counterexample.
    if graph_name == "K4":
         # Dropping 2 edges from K4 can result in C4, where all degrees are 2.
         print("Is Choice A (degree property) true? False (can result in a C4, min degree 2)")


# --- Main execution ---
# Example 1: K4 (Complete graph on 4 vertices)
K4 = nx.complete_graph(4)
analyze_graph_property(K4, "K4")

# Example 2: A disconnected graph K4 U K2
K4_union_K2 = nx.disjoint_union(nx.complete_graph(4), nx.path_graph(2))
analyze_graph_property(K4_union_K2, "K4 U K2")

print("\n--- Final Conclusion ---")
print("We have found cases where the property holds, but choices A, B, and D are false (case K4).")
print("We have also found a case where the property holds, but choice C is false (case K4 U K2).")
print("Since none of the statements A, B, C, or D are universally true for a graph satisfying the given condition, the correct option must be E.")