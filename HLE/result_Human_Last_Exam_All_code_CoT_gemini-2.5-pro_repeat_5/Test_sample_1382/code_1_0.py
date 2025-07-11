import networkx as nx
import numpy as np

def demonstrate_graph_property():
    """
    This function demonstrates the relationship c = lambda_n(G) / 2
    for a specific graph.
    """
    # 1. Create a graph that satisfies the property.
    # We'll use a graph with two disjoint components: a 4-cycle and a single edge (K2).
    # This graph has n=6 nodes and c=2 connected components.
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)]) # Component 1: 4-cycle
    G.add_edge(4, 5) # Component 2: K2 (a single edge)

    # 2. Calculate the number of connected components, c.
    num_components = nx.number_connected_components(G)
    
    # 3. Calculate the largest eigenvalue of the Laplacian, lambda_n.
    # The Laplacian L = D - A, where D is the degree matrix and A is the adjacency matrix.
    L = nx.laplacian_matrix(G).toarray()
    
    # The eigenvalues are calculated from the matrix L.
    # We need to use floating point numbers for precision.
    eigenvalues = np.linalg.eigvalsh(L.astype(float))
    lambda_n = max(eigenvalues)

    # 4. Verify the relationship from the problem.
    # The relationship is c = lambda_n / 2.
    # We check if the number of components we found is equal to half the largest eigenvalue.
    
    # The right hand side of the equation
    rhs = lambda_n / 2
    
    print("--- Verifying the Property ---")
    print(f"Graph: A disjoint union of a 4-cycle and a single edge.")
    print(f"Number of connected components (c): {num_components}")
    print(f"Largest eigenvalue of the Laplacian (lambda_n): {lambda_n:.4f}")
    print(f"Half of the largest eigenvalue (lambda_n / 2): {rhs:.4f}")
    
    print("\n--- Conclusion ---")
    print("The problem states that a certain property of the graph equals lambda_n / 2.")
    print("The most plausible interpretation is that this property is the number of connected components (c).")
    print("The statement tells us that the following equation holds for the graph:")
    print(f"Number of connected components = (Largest Laplacian Eigenvalue) / 2")
    print(f"In our example: {num_components} = {lambda_n:.2f} / 2")
    print("This corresponds to answer choice D.")


demonstrate_graph_property()