import networkx as nx

def find_maximal_induced_matching(graph):
    """
    Finds a maximal induced matching using a simple greedy algorithm.
    This is not guaranteed to be the maximum matching, but provides a lower bound.
    """
    induced_matching = []
    # We work on a copy of the graph to be able to remove nodes
    g = graph.copy()
    
    # While there are still edges in the graph
    while g.number_of_edges() > 0:
        # Pick an arbitrary edge
        u, v = list(g.edges())[0]
        induced_matching.append((u,v))
        
        # Get the set of vertices to remove: the endpoints and their neighbors
        nodes_to_remove = set([u, v])
        # Note: use original graph `graph` to get full neighborhood
        for node in [u, v]:
            nodes_to_remove.update(graph.neighbors(node))
            
        # Remove these nodes from our working copy
        g.remove_nodes_from(list(nodes_to_remove))
        
    return induced_matching

def main():
    """
    Main function to run the demonstration.
    """
    k = 20  # Dimension for the grid
    G = nx.grid_2d_graph(k, k)
    
    n = G.number_of_nodes()
    # For a grid graph, max degree is 4 (for non-boundary nodes)
    d = max(dict(G.degree()).values())

    print(f"Generated a {k}x{k} grid graph.")
    print(f"Number of vertices (n): {n}")
    print(f"Maximum degree (d): {d}")
    print("-" * 30)

    # Theoretical lower bound on induced matching size
    # From the bound alpha(G^2) >= n / (Delta(G^2) + 1) where Delta(G^2) <= d^2
    theoretical_bound = n / (d**2 + 1)
    
    # Find a maximal induced matching
    matching = find_maximal_induced_matching(G)
    matching_size = len(matching)

    print("Equation for theoretical lower bound of the induced matching size:")
    print(f"size >= n / (d^2 + 1)")
    print(f"size >= {n} / ({d}^2 + 1)")
    print(f"size >= {n} / ({d**2 + 1})")
    print(f"Theoretical lower bound: {theoretical_bound:.2f}")
    print("-" * 30)
    
    print(f"Found a maximal induced matching using a greedy algorithm.")
    print(f"Size of the found matching: {matching_size}")

    if matching_size >= theoretical_bound:
        print("\nAs shown, the size of the found matching respects the theoretical lower bound.")
        print("Since n can be arbitrarily large for a class C with unbounded treewidth,")
        print("the induced matching size must also be unbounded.")
    
if __name__ == '__main__':
    main()
