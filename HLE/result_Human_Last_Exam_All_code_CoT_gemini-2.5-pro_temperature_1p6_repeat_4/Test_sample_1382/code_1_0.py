import networkx as nx
import numpy as np

def analyze_graph_property(G):
    """
    Analyzes a graph to check the relationship proposed in answer choice D.

    Args:
        G (nx.Graph): A networkx graph object.
    """
    if G.number_of_nodes() <= 3:
        print("The graph must have more than 3 nodes.")
        return

    # Calculate the number of connected components, c
    c = nx.number_connected_components(G)

    # Calculate the largest eigenvalue of the graph Laplacian, lambda_n
    # Note: networkx returns eigenvalues in increasing order.
    # For a graph with c components, there will be c eigenvalues equal to 0.
    try:
        lambda_all = nx.laplacian_spectrum(G)
        lambda_n = lambda_all[-1]
    except Exception as e:
        print(f"Could not compute eigenvalues: {e}")
        return

    # The equation from choice D is c = lambda_n / 2
    # We will print the values of the left side (LHS) and right side (RHS)
    lhs = c
    rhs = lambda_n / 2

    print(f"For the given graph ({G.name}):")
    print(f"Number of nodes (n): {G.number_of_nodes()}")
    print(f"Number of edges (m): {G.number_of_edges()}")
    print("-" * 20)
    print("Checking the equation from Choice D: c = lambda_n(G) / 2")
    print(f"Number of connected components (c): {lhs}")
    print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")
    print(f"So, lambda_n(G) / 2 = {rhs:.4f}")
    print("-" * 20)
    print("The final equation is:")
    print(f"{lhs} = {rhs:.4f}")


# Example 1: A cycle graph on 5 nodes (C5).
# Here n=5, m=5, c=1. It has one component which is not bipartite.
G1 = nx.cycle_graph(5)
G1.name = "Cycle graph C5"
analyze_graph_property(G1)

print("\n" + "="*40 + "\n")

# Example 2: A disjoint union of two 4-cycles (C4 U C4).
# Here n=8, m=8, c=2. It has two components, both are bipartite.
G2 = nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
G2.name = "Disjoint union of two C4"
analyze_graph_property(G2)