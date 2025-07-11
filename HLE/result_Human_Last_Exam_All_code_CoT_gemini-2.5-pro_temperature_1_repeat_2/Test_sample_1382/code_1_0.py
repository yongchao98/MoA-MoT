import numpy as np
import networkx as nx

def check_k4_properties():
    """
    This function analyzes the K4 graph to verify the properties discussed.
    """
    # Create K4 graph
    G = nx.complete_graph(4)
    n = G.number_of_nodes()
    m = G.number_of_edges()

    # Get the unoriented incidence matrix B
    # Note: networkx incidence_matrix is oriented by default. We take the absolute value
    # to get the unoriented version as defined in the problem.
    B = np.abs(nx.incidence_matrix(G, oriented=True).toarray())

    # 1. Calculate null(B^T B) = null(B)
    # null(B) = m - rank(B)
    rank_B = np.linalg.matrix_rank(B)
    nullity_B = m - rank_B

    # 2. Calculate lambda_n(G) / 2
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = max(eigenvalues)
    
    # 3. Check the condition from the problem
    condition_holds = np.isclose(nullity_B, lambda_n / 2)

    print(f"Analyzing the graph K4 (complete graph on 4 nodes):")
    print(f"Number of nodes n = {n}")
    print(f"Number of edges m = {m}")
    print("-" * 30)

    # Check bipartiteness to find c_b
    # A graph is bipartite if it has no odd-length cycles.
    # K4 has C3, so it's not bipartite. Since it's connected, c_b = 0.
    c = nx.number_connected_components(G)
    c_b = 0 # K4 is connected and not bipartite
    
    print(f"Number of connected components c = {c}")
    print(f"Number of bipartite components c_b = {c_b}")
    
    # Verify our formula for nullity
    # Formula: nullity = m - n + c_b
    formula_nullity = m - n + c_b
    print(f"Nullity of B calculated from rank: null(B) = {nullity_B}")
    print(f"Nullity of B calculated from formula (m-n+c_b): {formula_nullity}")
    print("-" * 30)
    
    # Verify the problem's core statement
    print(f"The problem states: null(B^T B) = lambda_n(G) / 2")
    print(f"For K4, null(B^T B) = {nullity_B}")
    print(f"For K4, lambda_n(G) = {lambda_n:.4f}")
    print(f"So, lambda_n(G) / 2 = {lambda_n / 2:.4f}")
    print(f"Does the condition hold? {condition_holds}")
    print("-" * 30)

    # 4. Check the options
    print("Evaluating the answer choices for K4:")
    
    # Option A
    # Drop lambda_n/2 = 2 edges. Does it result in >= 2 nodes with degree <= 1?
    # Remove two non-adjacent edges (0,1) and (2,3) to create a C4
    G_A = G.copy()
    G_A.remove_edges_from([(0,1), (2,3)])
    degrees_A = [d for n, d in G_A.degree()]
    nodes_le_1 = sum(1 for d in degrees_A if d <= 1)
    print(f"Option A Check: Drop {lambda_n/2} edges.")
    print(f"  Removing edges (0,1) and (2,3) results in degrees: {degrees_A}")
    print(f"  Number of nodes with degree <= 1 is {nodes_le_1}.")
    print(f"  Is {nodes_le_1} >= 2? {nodes_le_1 >= 2}. So Option A is FALSE.")

    # Option B
    print(f"Option B Check: Is c >= lambda_n/2?")
    print(f"  Is {c} >= {lambda_n/2}? {c >= lambda_n/2}. So Option B is FALSE.")

    # Option D
    print(f"Option D Check: Is c == lambda_n/2?")
    print(f"  Is {c} == {lambda_n/2}? {c == lambda_n/2}. So Option D is FALSE.")

check_k4_properties()