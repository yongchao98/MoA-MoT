import collections

def demonstrate_non_uniqueness():
    """
    This function demonstrates that the 1-crossing drawing is not unique by
    modeling a specific maximal planar graph and showing the options for crossing.
    """
    # Define a maximal planar graph G on 5 vertices {0, 1, 2, 3, 4}.
    # Construction: A central triangle (1,2,3) separates an outer vertex 0
    # and an inner vertex 4.
    # Both 0 and 4 are connected to all of {1, 2, 3}.
    G = collections.defaultdict(set)
    edges = [
        (0, 1), (0, 2), (0, 3),  # Outer vertex connections
        (4, 1), (4, 2), (4, 3),  # Inner vertex connections
        (1, 2), (2, 3), (3, 1)   # The separating triangle
    ]
    for u, v in edges:
        G[u].add(v)
        G[v].add(u)
    
    # The edge to add, 'e', connects the only two non-adjacent vertices.
    e = (0, 4)
    u, v = e
    
    print(f"Consider a maximal planar graph G on vertices {sorted(list(G.keys()))}.")
    print(f"Let the new edge to be added be e = {e}.")
    print("\nIn a planar embedding of G, vertices 1, 2, and 3 form a cycle that")
    print(f"separates vertex {u} from vertex {v}.")
    
    # The edges of the separating cycle
    separating_cycle_edges = [(1, 2), (2, 3), (3, 1)]
    
    print("\nTo draw the new edge e = (0, 4), we must cross an edge of this separating cycle.")
    print("This can be done with a single crossing. However, the choice of which edge to cross is not unique.")
    print("\nHere are the possible pairs of edges that can form the single crossing:")
    
    # Print the different options for the crossing
    for i, crossed_edge in enumerate(separating_cycle_edges):
        # We print each number in the final equation as requested.
        print(f"Option {i+1}: New edge ({e[0]}, {e[1]}) crosses existing edge ({crossed_edge[0]}, {crossed_edge[1]})")

    print("\nSince there are multiple distinct choices for the single crossing,")
    print("the drawing of G' = G U {e} with at most one crossing is not unique.")

demonstrate_non_uniqueness()
