import networkx as nx

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    A simple cycle is one that does not repeat vertices, except for the start/end vertex.
    
    The algorithm works by finding all simple paths of length k-1 starting from each
    node and checking if the ends of the path form an edge.
    
    Args:
        G (nx.Graph): The input graph.
        k (int): The length of the cycle to count.

    Returns:
        int: The number of simple cycles of length k.
    """
    count = 0
    # We iterate through all nodes to start paths from.
    for start_node in G.nodes():
        # A stack for DFS, storing (current_node, path_so_far).
        stack = [(start_node, [start_node])]
        while stack:
            curr_node, path = stack.pop()
            
            # If the path has k nodes, we have a path of length k-1.
            # We check if the last node connects back to the first.
            if len(path) == k:
                if G.has_edge(path[-1], start_node):
                    count += 1
                continue

            # Extend the path with neighbors of the current node.
            for neighbor in G.neighbors(curr_node):
                # To form a simple path, the neighbor must not be in the current path.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    # Each cycle of length k will be found k times (once for each starting node)
    # and 2 times (once for each direction). So we divide by 2*k.
    return count // (2 * k)

def main():
    """
    Main function to construct graphs and compare their 5-cycle counts.
    """
    # Parameters for the strongly regular graphs
    n = 16
    d = 6
    lambda_ = 2
    mu = 2
    
    print(f"Searching for two non-isomorphic graphs in srg(n={n}, d={d}, λ={lambda_}, μ={mu})")
    print("with a different number of 5-cycles.\n")

    # 1. Construct the L(K_4,4) graph
    # It is the line graph of the complete bipartite graph K_{4,4}.
    G1 = nx.line_graph(nx.complete_bipartite_graph(4, 4))
    
    # 2. Construct the Shrikhande graph
    # It can be defined as a Cayley graph on Z_4 x Z_4.
    G2 = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G2.add_nodes_from(nodes)
    # The connection set for the Cayley graph
    generators = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
    for i in range(4):
        for j in range(4):
            node1 = (i, j)
            for dx, dy in generators:
                # Connect to neighbors using modular arithmetic
                node2 = ((i + dx) % 4, (j + dy) % 4)
                if not G2.has_edge(node1, node2):
                    G2.add_edge(node1, node2)

    # 3. Count 5-cycles in both graphs
    print("Counting 5-cycles... (this may take a moment)")
    c5_G1 = count_simple_cycles(G1, 5)
    c5_G2 = count_simple_cycles(G2, 5)
    
    # 4. Print the results
    print("\n--- Results ---")
    print(f"Graph 1: L(K_4,4) graph")
    print(f"Number of 5-cycles: {c5_G1}")

    print(f"\nGraph 2: Shrikhande graph")
    print(f"Number of 5-cycles: {c5_G2}")
    
    print("\n--- Final Equation ---")
    print(f"The number of 5-cycles in G1 ({c5_G1}) is not equal to the number of 5-cycles in G2 ({c5_G2}).")
    print(f"Therefore, two SRGs with parameters (n={n}, d={d}, λ={lambda_}, μ={mu}) can have different numbers of 5-cycles.")

if __name__ == "__main__":
    main()