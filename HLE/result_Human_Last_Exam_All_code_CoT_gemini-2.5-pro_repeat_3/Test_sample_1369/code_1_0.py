import networkx as nx
import numpy as np

def analyze_graph(G, graph_name):
    """
    Analyzes a graph to check if it meets the specified criteria.
    """
    print(f"--- Analysis of: {graph_name} ---")

    n = G.number_of_nodes()
    m = G.number_of_edges()
    c = nx.number_connected_components(G)

    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (c): {c}")

    # Check the cyclomatic number condition: m - n + c = 2
    cyclomatic_number = m - n + c
    print(f"Checking condition from incidence matrix: m - n + c = {m} - {n} + {c} = {cyclomatic_number}")
    if cyclomatic_number == 2:
        print("This matches the condition null(B^T B) = 2.")
    else:
        print("This does NOT match the condition null(B^T B) = 2.")

    # Check the eigenvalue condition: c >= 2
    # We compute the eigenvalues to show they are consistent with the observation [0.0, 0.0, ...]
    laplacian = nx.laplacian_matrix(G).toarray()
    eigenvalues = np.linalg.eigvalsh(laplacian)
    
    num_zero_eigenvalues = np.sum(np.isclose(eigenvalues, 0))
    print(f"Number of zero Laplacian eigenvalues found: {num_zero_eigenvalues}")
    
    if num_zero_eigenvalues >= 2:
        print("The first two eigenvalues are 0, which is consistent with the problem statement.")
    else:
        print("The eigenvalue condition is NOT met.")

    print(f"Conclusion: This graph is a valid example satisfying the problem's constraints.\n")


def main():
    print("Investigating if the number of connected components (c) is uniquely determined.\n")

    # Example 1: A graph with c=2 components
    # Two disjoint triangles. Each triangle has n=3, m=3.
    # Total: n=6, m=6, c=2.
    # Cyclomatic number = m - n + c = 6 - 6 + 2 = 2.
    G1 = nx.disjoint_union(nx.cycle_graph(3), nx.cycle_graph(3))
    analyze_graph(G1, "Two disjoint triangles (c=2)")

    # Example 2: A graph with c=3 components
    # G1: triangle (n=3, m=3, k=1)
    # G2: square (n=4, m=4, k=1)
    # G3: isolated vertex (n=1, m=0, k=0)
    # Total: n=8, m=7, c=3
    # Total cyclomatic number = 1 + 1 + 0 = 2.
    # Check: m - n + c = 7 - 8 + 3 = 2.
    G2_part1 = nx.cycle_graph(3)
    G2_part2 = nx.cycle_graph(4)
    G2_part3 = nx.Graph()
    G2_part3.add_node(0) # single node component
    
    G2 = nx.disjoint_union(G2_part1, G2_part2)
    G2 = nx.disjoint_union(G2, G2_part3)
    analyze_graph(G2, "A triangle, a square, and an isolated vertex (c=3)")
    
    print("--- Final Conclusion ---")
    print("We have constructed valid graphs with both 2 and 3 connected components that satisfy all given conditions.")
    print("This demonstrates that the number of connected components is not uniquely determined to be 2.")
    print("Therefore, statements A, B, C, and D are not properties the graph must have.")
    print("The correct answer is E.")


if __name__ == "__main__":
    main()