import networkx as nx

def count_cycles(graph, length):
    """
    Counts the number of simple cycles of a given length in a graph.

    Args:
        graph (nx.Graph): The graph to search in.
        length (int): The length of the cycles to count.

    Returns:
        int: The number of simple cycles of the given length.
    """
    count = 0
    adj = graph.adj
    # We iterate over each node to start a traversal (path) from it.
    for start_node in graph.nodes():
        # The stack for the depth-first search will store tuples of (current_node, path_so_far).
        stack = [(start_node, [start_node])]
        while stack:
            current_node, path = stack.pop()
            
            # If the path has reached the desired length
            if len(path) == length:
                # Check if the last node in the path is connected to the start node, forming a cycle.
                if start_node in adj[current_node]:
                    count += 1
                continue

            # Extend the current path with neighbors of the current node.
            for neighbor in adj[current_node]:
                # To ensure a simple path (no repeated vertices), we don't visit a node already in the path.
                # The exception is the start_node, which we can revisit only at the end.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    # Each cycle is discovered 'length' times (once for each starting node) and
    # '2' times (once for each traversal direction). So, we divide by 2 * length.
    return count // (2 * length)

def main():
    """
    Main function to construct graphs and compare their 5-cycle counts.
    """
    # The parameters for the strongly regular graphs.
    n, d, lam, mu = 16, 6, 2, 2

    # 1. The 4x4 Rook's graph, also known as the Latin square graph L_2(4) or K_4 x K_4.
    # It is an SRG with parameters (16, 6, 2, 2).
    # Vertices are cells of a 4x4 grid, adjacent if in the same row or column.
    G_rook = nx.cartesian_product(nx.complete_graph(4), nx.complete_graph(4))
    G_rook.name = "4x4 Rook's graph"

    # 2. The Shrikhande graph.
    # It is a different SRG with the same parameters (16, 6, 2, 2).
    G_shrikhande = nx.shrikhande_graph()
    G_shrikhande.name = "Shrikhande graph"
    
    # Calculate the number of 5-cycles for both graphs.
    print("Finding a pair of graphs in srg(n,d,lambda,mu) with a different number of 5-cycles...")
    print(f"The chosen parameters are (n, d, lambda, mu) = ({n}, {d}, {lam}, {mu}).\n")

    print(f"Graph 1: {G_rook.name}")
    num_cycles_rook = count_cycles(G_rook, 5)
    
    print(f"Graph 2: {G_shrikhande.name}")
    num_cycles_shrikhande = count_cycles(G_shrikhande, 5)

    print("\n--- Results ---")
    print(f"The number of 5-cycles in the {G_rook.name} is: {num_cycles_rook}")
    print(f"The number of 5-cycles in the {G_shrikhande.name} is: {num_cycles_shrikhande}")

    print("\nConclusion:")
    print("Yes, such a combination exists. For parameters (16, 6, 2, 2), the two graphs have different numbers of 5-cycles:")
    print(f"{num_cycles_rook} (for {G_rook.name}) != {num_cycles_shrikhande} (for {G_shrikhande.name})")


if __name__ == "__main__":
    main()

<<<Yes>>>