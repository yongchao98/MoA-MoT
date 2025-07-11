import numpy as np
import networkx as nx

def analyze_graph_property(G, graph_name):
    """
    Analyzes a graph based on the properties mentioned in the problem.
    """
    print(f"--- Analyzing graph: {graph_name} ---")

    # n = number of nodes, m = number of edges
    n = G.number_of_nodes()
    m = G.number_of_edges()

    # k = number of connected components
    k = nx.number_connected_components(G)

    # L = Laplacian matrix
    L = nx.laplacian_matrix(G).toarray()
    
    # Calculate eigenvalues of L and find the largest, lambda_n
    # We use eigvalsh because the Laplacian is symmetric, which is more stable
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    # Calculate k_bip = number of bipartite connected components
    k_bip = 0
    for component_nodes in nx.connected_components(G):
        subgraph = G.subgraph(component_nodes)
        if nx.is_bipartite(subgraph):
            k_bip += 1

    # Calculate the quantities from the problem
    # LHS of the given equation: null(B^T B) = m - n + k_bip
    # We will call the number of components k = null(L)
    lhs_original = m - n + k_bip
    rhs = lambda_n / 2
    
    # The number of connected components, k
    num_components = k
    
    # ---- The Reasoning ----
    # The problem gives the equation: null(B^T B) = lambda_n / 2
    # The term null(B^T B) for an unoriented incidence matrix is m - n + k_bip.
    # The answers, however, relate to the number of connected components, k.
    # The most fundamental result in spectral graph theory is that the number of
    # connected components 'k' is the nullity of the Laplacian matrix, L.
    # null(L) = k
    # It is highly likely that the convoluted phrase "null(B^T B)" was a mistake
    # for the much more fundamental quantity representing components, which is 'k'.
    # So, "what they are trying to tell you" is that for this graph, k = lambda_n / 2.
    # This statement is precisely Answer D.

    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (k): {k}")
    print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")

    print("\nInterpreting the statement 'k = lambda_n / 2':")
    # This represents the implication of the statement, which is answer D
    # We will output the numbers that form the equation of Answer D
    print(f"The number of components is {num_components}.")
    print(f"The value of lambda_n / 2 is {rhs:.4f}.")
    print("Answer D suggests the final equation:")
    print(f"{num_components} = {lambda_n:.4f} / 2")
    print("\n")


# Let's create an example graph that actually satisfies the original obscure condition
# This graph consists of two disjoint squares (4-cycles)
G1 = nx.cycle_graph(4)
G2 = nx.cycle_graph(4)
G = nx.disjoint_union(G1, G2)
analyze_graph_property(G, "Two Disjoint Squares (C4 U C4)")

# Another example: A single 4-cycle graph
G_C4 = nx.cycle_graph(4)
analyze_graph_property(G_C4, "Single Square (C4)")


<<<D>>>