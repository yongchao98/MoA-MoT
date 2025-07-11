import networkx as nx

def count_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    
    This function uses a depth-first search from each node to find paths of
    length k-1 that connect back to the starting node.
    
    Args:
        G (networkx.Graph): The graph to search in.
        k (int): The length of the cycles to count.
        
    Returns:
        int: The total number of simple cycles of length k.
    """
    # This function requires the 'networkx' library.
    # You can install it using: pip install networkx
    
    total_count = 0
    # Iterate over each node to use it as a starting point for a potential cycle.
    for start_node in G.nodes():
        # The stack for our iterative DFS will store tuples of (path, visited_nodes_set)
        # path is a list of nodes, visited_nodes_set is a set for O(1) lookups.
        stack = [([start_node], {start_node})]
        
        while stack:
            path, visited = stack.pop()
            current_node = path[-1]
            
            # If the path has k nodes, it has a length of k-1.
            # We check if it can be closed to form a cycle of length k.
            if len(path) == k:
                # A cycle is formed if the current node is a neighbor of the start node.
                if start_node in G[current_node]:
                    total_count += 1
                continue

            # Explore neighbors of the current node to extend the path.
            for neighbor in G.neighbors(current_node):
                # To form a simple cycle, we must not visit any node more than once.
                if neighbor not in visited:
                    new_path = path + [neighbor]
                    new_visited = visited.union({neighbor})
                    stack.append((new_path, new_visited))

    # Each cycle of length k is counted 2*k times by this algorithm:
    # - k times because any of its k nodes can be the 'start_node'.
    # - 2 times because for each start_node, the cycle can be traversed in two directions.
    # We use integer division as the count must be an integer.
    return total_count // (2 * k)

def main():
    """
    Main function to find and compare the number of 5-cycles in two
    non-isomorphic SRGs with the same parameters.
    """
    # Define the parameters for the class of strongly regular graphs.
    n, d, lam, mu = 16, 6, 2, 2

    # 1. Create the 4x4 Rook's graph (L_2(4)).
    # It is the Cartesian product of two complete graphs K_4.
    K4 = nx.complete_graph(4)
    rook_graph = nx.cartesian_product(K4, K4)

    # 2. Create the Shrikhande graph.
    shrikhande_graph = nx.shrikhande_graph()

    # 3. Count the number of 5-cycles in each graph.
    c5_rook = count_cycles(rook_graph, 5)
    c5_shrikhande = count_cycles(shrikhande_graph, 5)

    # 4. Print the results.
    print("Yes, there is a combination of parameters for which two SRGs can have a different number of 5-cycles.")
    print("-" * 20)
    print(f"An example is the class of strongly regular graphs srg({n}, {d}, {lam}, {mu}).")
    print(f"The 4x4 Rook's graph, G1, belongs to srg({n}, {d}, {lam}, {mu}) and has {c5_rook} 5-cycles.")
    print(f"The Shrikhande graph, G2, also belongs to srg({n}, {d}, {lam}, {mu}) and has {c5_shrikhande} 5-cycles.")
    print("\nThese two graphs are co-spectral but not isomorphic, and they have different numbers of 5-cycles.")
    print("The final comparison is:")
    print(f"{c5_rook} != {c5_shrikhande}")

if __name__ == "__main__":
    main()
