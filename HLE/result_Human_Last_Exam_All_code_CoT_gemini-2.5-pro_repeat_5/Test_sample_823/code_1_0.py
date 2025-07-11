import networkx as nx

def find_maximal_induced_matching(G):
    """
    Finds a maximal induced matching in a graph G using a greedy algorithm.
    An induced matching is a set of edges with no common vertices and no
    edges between the vertices of any two edges in the set.
    """
    # Work on a copy to not modify the original graph object
    g_copy = G.copy()
    matching = []
    
    # Iterate as long as there are edges in the graph
    while g_copy.number_of_edges() > 0:
        # Pick an arbitrary edge from the remaining graph
        u, v = next(iter(g_copy.edges()))
        
        # Add the chosen edge to our matching
        matching.append((u, v))
        
        # Identify all vertices to be removed.
        # This includes the edge's endpoints and all of their neighbors
        # in the *original* graph G.
        nodes_to_remove = {u, v}
        nodes_to_remove.update(G.neighbors(u))
        nodes_to_remove.update(G.neighbors(v))
        
        # Remove these vertices from our working copy of the graph
        g_copy.remove_nodes_from(nodes_to_remove)
        
    return matching

def demonstrate_induced_matching():
    """
    Demonstrates that a large graph with bounded degree has a large
    induced matching. We use a grid graph as an example.
    """
    k = 30  # Side length of the grid
    G = nx.grid_2d_graph(k, k)
    
    n = G.number_of_nodes()
    
    # Max degree of an infinite grid is 4. For a finite kxk grid, it's 4 (for k>2).
    # We can calculate it directly for confirmation.
    d = max(dict(G.degree()).values())

    print(f"Generated a {k}x{k} grid graph.")
    print(f"Number of vertices (n): {n}")
    print(f"Maximum degree (d): {d}")
    print("-" * 30)

    # Find a maximal induced matching
    maximal_matching = find_maximal_induced_matching(G)
    matching_size = len(maximal_matching)
    
    # Calculate the theoretical lower bound for the matching size
    lower_bound = n / (2 * (d + 1))
    
    print("Finding a maximal induced matching...")
    print(f"Found a matching of size: {matching_size}")
    
    print("\nChecking the theoretical bound:")
    print(f"|M| >= n / (2 * (d + 1))")
    print("Substituting the values:")
    print(f"{matching_size} >= {n} / (2 * ({d} + 1))")
    print(f"{matching_size} >= {n} / {2 * (d + 1)}")
    print(f"{matching_size} >= {lower_bound:.2f}")

    if matching_size >= lower_bound:
        print("\nThe found matching size satisfies the theoretical lower bound.")
    else:
        print("\nError: The found matching size is smaller than the theoretical bound.")

if __name__ == '__main__':
    demonstrate_induced_matching()
