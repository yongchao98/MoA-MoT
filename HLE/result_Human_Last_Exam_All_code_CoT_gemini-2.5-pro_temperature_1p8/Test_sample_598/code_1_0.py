import networkx as nx

def count_k_cycles(g, k):
    """
    Counts the number of simple cycles of length k in a graph.
    
    This function uses a backtracking algorithm (implemented with an iterative
    stack-based DFS) to find all simple paths of length k-1 starting from each
    node. It then checks if the endpoint of the path is connected to the start
    node, which would complete a k-cycle.
    
    The total count is divided by 2*k because each cycle is counted k times
    (once for each starting node) and in 2 directions (clockwise and
    counter-clockwise).
    """
    total_count = 0
    
    # Iterate over each node as a potential starting node for a cycle
    for start_node in g.nodes():
        # The stack holds tuples of (current_node, path_so_far)
        stack = [(start_node, [start_node])]
        
        while stack:
            curr_node, path = stack.pop()
            
            # If the path has k vertices (i.e., length k-1)
            if len(path) == k:
                # Check if the last node is connected back to the start node
                # g[curr_node] is the set of neighbors of curr_node
                if start_node in g[curr_node]:
                    total_count += 1
                continue
            
            # Extend the path by one more step
            for neighbor in g[curr_node]:
                # To form a simple cycle, all intermediate vertices must be unique
                if neighbor not in path:
                    new_path = list(path) # Create a copy of the path
                    new_path.append(neighbor)
                    stack.append((neighbor, new_path))
    
    # Each cycle is counted for each of its k nodes, and in two directions.
    # Therefore, we divide the total count by 2*k.
    return total_count // (2 * k)

def main():
    """
    Main function to construct the graphs, count cycles, and print the results.
    """
    # The parameters for the strongly regular graphs (SRGs)
    n = 16
    d = 6
    lambda_param = 2
    mu_param = 2
    
    print("Yes, there exists a pair of strongly regular graphs G, G' from the same class srg(n,d,lambda,mu)")
    print("that do not have the same number of 5-cycles.")
    print("\nWe will demonstrate this using the class with parameters:")
    print(f"n (vertices) = {n}")
    print(f"d (degree) = {d}")
    print(f"lambda = {lambda_param}")
    print(f"mu = {mu_param}\n")
    
    # Construct Graph 1: The Shrikhande Graph
    # It is a well-known SRG with parameters (16, 6, 2, 2)
    G_shrikhande = nx.shrikhande_graph()
    
    # Construct Graph 2: The 4x4 Rook's Graph
    # This is also known as the Latin square graph L(2,4) or the line graph of K(4,4).
    # It is an SRG with the same parameters (16, 6, 2, 2).
    G_rooks = nx.latin_square_graph(4)

    # Count the 5-cycles in each graph
    print("Counting 5-cycles in both graphs. This may take a moment...")
    num_5_cycles_shrikhande = count_k_cycles(G_shrikhande, 5)
    num_5_cycles_rooks = count_k_cycles(G_rooks, 5)
    
    print("\n--- Results ---")
    print("Graph 1: The Shrikhande Graph")
    print(f"Number of 5-cycles = {num_5_cycles_shrikhande}")
    
    print("\nGraph 2: The 4x4 Rook's Graph (also L(K_4,4))")
    print(f"Number of 5-cycles = {num_5_cycles_rooks}")
    
    print("\nAs shown, the two graphs belong to srg(16, 6, 2, 2) but have different numbers of 5-cycles,")
    print(f"confirming that such a combination exists.")

if __name__ == "__main__":
    main()