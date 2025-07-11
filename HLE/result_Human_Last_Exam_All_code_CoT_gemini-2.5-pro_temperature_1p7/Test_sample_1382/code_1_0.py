import networkx as nx
import numpy as np

def analyze_graph_property(G):
    """
    Analyzes a graph based on the property null(B^T B) = lambda_n(G) / 2.
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    
    if n <= 3:
        print("Graph must have more than 3 nodes.")
        return

    print(f"Graph properties: Nodes n={n}, Edges m={m}")

    # For the logic of option A, we need the number of connected components.
    # networkx.number_connected_components is efficient.
    k = nx.number_connected_components(G)
    print(f"Number of connected components k={k}")

    # --- Step 1: Verify the condition ---

    # Calculate the left side: null(B^T B)
    # This is equivalent to the cyclomatic number m - n + k.
    cyclomatic_number = m - n + k
    
    # Calculate the right side: lambda_n(G) / 2
    # Use the normalized laplacian to avoid scaling issues if desired, but standard works fine.
    # networkx returns a Scipy sparse matrix, so convert to dense array for numpy's eigvalsh
    L = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1] # Eigenvalues from eigvalsh are sorted
    
    # Using np.isclose for robust floating point comparison
    condition_holds = np.isclose(cyclomatic_number, lambda_n / 2)
    
    print(f"\nVerifying the condition: null(B^T B) = lambda_n(G) / 2")
    print(f"Left side (cyclomatic number m-n+k): {m} - {n} + {k} = {cyclomatic_number}")
    print(f"Right side (lambda_n / 2): {lambda_n:.4f} / 2 = {lambda_n/2:.4f}")
    
    if not condition_holds:
        print("The given graph does not satisfy the condition.")
        return
    else:
        print("The condition holds for this graph.")

    # --- Step 2: Test the implication from Choice A ---
    
    num_edges_to_drop = int(round(lambda_n / 2))
    print(f"\nTesting choice A: Drop {num_edges_to_drop} edges to form a forest.")

    # A graph's cycle basis has a length equal to the cyclomatic number.
    # Removing one edge from each cycle in the basis makes the graph acyclic (a forest).
    G_prime = G.copy()
    cycles = nx.cycle_basis(G_prime)
    
    # Note: cycle_basis finds fundamental cycles. We need to remove 'cyclomatic_number' edges.
    edges_dropped = 0
    if not cycles:
        print("The graph is already a forest.")
    else:
        # Remove one edge from each fundamental cycle found
        for cycle in cycles:
            if edges_dropped < num_edges_to_drop:
                 # Remove the first edge of the cycle
                G_prime.remove_edge(cycle[0], cycle[1])
                edges_dropped += 1

    print(f"Created a new graph G' by dropping {edges_dropped} edge(s).")
    print(f"Is G' a forest? {nx.is_forest(G_prime)}")
    
    # Check the degrees in the new graph G'
    degrees = [d for n, d in G_prime.degree()]
    nodes_with_low_degree = [node for node, degree in G_prime.degree() if degree <= 1]
    
    print(f"Degrees in G': {degrees}")
    print(f"Number of nodes with degree <= 1 in G': {len(nodes_with_low_degree)}")
    
    if len(nodes_with_low_degree) >= 2:
        print("Conclusion of choice A holds true: The resulting forest has at least two nodes with degree <= 1.")
    else:
        print("Conclusion of choice A does not hold.")


# Create the counterexample graph G = C_4 U C_4 which is disconnected but satisfies the property.
C4_1 = nx.cycle_graph(4)
C4_2 = nx.cycle_graph(4)
# Relabel nodes of the second C4 to make them distinct for the disjoint union
C4_2 = nx.relabel_nodes(C4_2, {0:4, 1:5, 2:6, 3:7})
G = nx.disjoint_union(C4_1, C4_2)

analyze_graph_property(G)
