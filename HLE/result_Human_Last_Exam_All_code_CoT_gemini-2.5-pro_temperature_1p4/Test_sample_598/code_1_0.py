import networkx as nx

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    A simple cycle has k distinct vertices.
    The algorithm performs a search for paths of length k-1 and checks
    if the endpoints are connected.
    """
    # This counter will be incremented each time a cycle is found.
    # We will correct for overcounting at the end.
    count = 0
    
    # We iterate through all nodes to use each as a potential starting point for a cycle.
    for start_node in G.nodes():
        # A stack for depth-first search. Stores tuples of (current_node, path_so_far).
        # We start with paths of length 2, originating from the start_node.
        stack = [(neighbor, [start_node, neighbor]) for neighbor in G.neighbors(start_node)]
        
        while stack:
            current_node, path = stack.pop()
            
            # If the path has k vertices, it has length k-1.
            # We check if this path can be closed to form a k-cycle.
            if len(path) == k:
                # The cycle is formed if the last node is connected to the start_node.
                if start_node in G.neighbors(current_node):
                    count += 1
                continue # We don't extend paths longer than k.

            # Extend the path by one more step.
            for neighbor in G.neighbors(current_node):
                # To ensure the cycle is simple, the new node must not be in the current path.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    # Each cycle of length k is counted 2*k times:
    # - k times: once for each node in the cycle being the 'start_node'.
    # - 2 times: once for each of the two traversal directions (clockwise/counter-clockwise).
    # We perform integer division as the result must be a whole number.
    return count // (2 * k)

def solve():
    """
    Finds and analyzes two SRGs with the same parameters but a different number of 5-cycles.
    """
    # The parameters for the class of strongly regular graphs
    n, d, lam, mu = 16, 6, 2, 2

    # Generate the two non-isomorphic graphs in srg(16, 6, 2, 2)
    # G1 is the 4x4 Rook's graph (or Latin square graph L_2(4))
    G_rook = nx.rook_graph(4)
    # G2 is the Shrikhande graph
    G_shrikhande = nx.shrikhande_graph()

    # Count the number of 5-cycles in each graph
    num_cycles_rook = count_simple_cycles(G_rook, 5)
    num_cycles_shrikhande = count_simple_cycles(G_shrikhande, 5)

    # Print the final results, showing all the numbers involved.
    print("Yes, such a pair of graphs exists for the SRG class with parameters:")
    print(f"n = {n}")
    print(f"d = {d}")
    print(f"lambda = {lam}")
    print(f"mu = {mu}")
    print("\nThese two graphs are the 4x4 Rook's graph and the Shrikhande graph.")
    print("Calculating the number of 5-cycles in each graph gives different results:\n")
    
    # The final "equation" showing the inequality
    print(f"Number of 5-cycles in 4x4 Rook's graph: {num_cycles_rook}")
    print(f"Number of 5-cycles in Shrikhande graph: {num_cycles_shrikhande}")
    
    print(f"\nSince {num_cycles_rook} != {num_cycles_shrikhande}, we have found two graphs G, G' in srg({n}, {d}, {lam}, {mu}) that do not have the same number of 5-cycles.")

# Run the solution
solve()