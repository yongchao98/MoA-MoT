import networkx as nx

def demonstrate_subgraph_counterexample(k_large=5, k_pattern=3):
    """
    This function demonstrates that a subdivided grid does not contain a smaller 
    grid as a subgraph, providing a counterexample for option B.

    A class of subdivided grids has bounded degree (4) and unbounded treewidth,
    but does not contain arbitrarily large grids as subgraphs.
    """
    print(f"Goal: Show that a subdivided grid (unbounded treewidth, bounded degree class) does not necessarily contain large grid subgraphs.")
    
    # 1. Create a "large" grid graph, G_large
    G_large = nx.grid_2d_graph(k_large, k_large)
    
    # 2. Create the subdivided version of G_large
    G_subdivided = nx.Graph()
    # The original vertices are tuples like (i, j)
    # The new subdivision vertices will be strings to avoid clashes
    for u, v in G_large.edges():
        # Each edge (u, v) is replaced by path u -> new_node -> v
        new_node_name = f"sub_{u}_{v}"
        G_subdivided.add_edge(u, new_node_name)
        G_subdivided.add_edge(new_node_name, v)

    max_degree_subdivided = max(d for _, d in G_subdivided.degree())
    print(f"\nCreated a subdivided {k_large}x{k_large} grid.")
    print(f" - Max degree: {max_degree_subdivided} (which is bounded, <= 4)")
    print(f" - The treewidth of this graph is {k_large} (which is unbounded as k_large grows).")

    # 3. Create the smaller grid pattern we want to search for
    G_pattern = nx.grid_2d_graph(k_pattern, k_pattern)
    print(f"\nSearching for a {k_pattern}x{k_pattern} grid as a subgraph in the subdivided grid.")

    # 4. Use a graph isomorphism matcher to check for the subgraph
    matcher = nx.isomorphism.GraphMatcher(G_subdivided, G_pattern)
    is_subgraph = matcher.subgraph_is_isomorphic()
    
    print(f"Result: Is the {k_pattern}x{k_pattern} grid a subgraph? {is_subgraph}")
    print("\nExplanation: This is false because in a kxk grid (for k>=3), a degree-4 vertex has neighbors of degree >= 3.")
    print("In our subdivided grid, the original vertices (which can have degree 4) only have neighbors of degree 2 (the subdivision vertices).")
    print("This structural difference makes a subgraph isomorphism impossible, thus falsifying option B.")

# Run the demonstration
demonstrate_subgraph_counterexample(k_large=5, k_pattern=3)