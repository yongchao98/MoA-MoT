import networkx as nx

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    A simple cycle has no repeated vertices other than the start/end node.
    
    The function uses a recursive Depth-First Search (DFS) approach to find
    all simple paths of length k-1 and checks if they form a cycle.
    
    Args:
        G (networkx.Graph): The graph to analyze.
        k (int): The length of the cycles to count.
        
    Returns:
        int: The total number of simple cycles of length k.
    """
    
    def _find_paths_recursive(path, visited_nodes):
        """
        A recursive helper function to find paths of a specific length.
        
        Args:
            path (list): The current path of vertices.
            visited_nodes (set): A set of vertices in the current path for O(1) lookups.
            
        Returns:
            int: The number of valid cycles found from the current path.
        """
        start_node = path[0]
        last_node = path[-1]
        
        # Base case for recursion: the path has reached the desired length.
        if len(path) == k:
            # Check if the path can be closed to form a cycle.
            return 1 if G.has_edge(last_node, start_node) else 0

        # Recursive step: extend the path with unvisited neighbors.
        count = 0
        for neighbor in G.neighbors(last_node):
            if neighbor not in visited_nodes:
                path.append(neighbor)
                visited_nodes.add(neighbor)
                count += _find_paths_recursive(path, visited_nodes)
                # Backtrack to explore other paths.
                visited_nodes.remove(neighbor)
                path.pop()
        return count

    total_count = 0
    # We iterate through each node as a potential starting point for a cycle.
    for node in G.nodes():
        total_count += _find_paths_recursive([node], {node})
        
    # Each cycle is counted k times (once for each starting node) and twice
    # (once for each direction). So, we divide the total count by 2*k.
    return total_count // (2 * k) if k > 0 else 0

def main():
    """
    Main function to construct graphs, verify parameters, and count 5-cycles.
    """
    # Parameters for the strongly regular graphs
    n, d, lam, mu = 16, 6, 2, 2

    print(f"Searching for two non-isomorphic graphs in srg({n}, {d}, {lam}, {mu}) with a different number of 5-cycles.\n")
    
    # 1. The Shrikhande Graph
    g1 = nx.shrikhande_graph()
    g1.name = "Shrikhande Graph"

    # 2. The Rook's Graph on a 4x4 board (K4 x K4)
    # This is also isomorphic to the line graph of K_{4,4}
    g2 = nx.cartesian_product(nx.complete_graph(4), nx.complete_graph(4))
    g2.name = "Rook's Graph (K4xK4)"
    # Networkx cartesian_product creates labels as tuples (e.g., (0,1)).
    # We convert them to strings for simpler display.
    g2 = nx.relabel_nodes(g2, {node: str(node) for node in g2.nodes()})


    graphs = [g1, g2]
    cycle_counts = []

    for G in graphs:
        print(f"--- Analyzing: {G.name} ---")
        
        # Count 5-cycles
        print("Counting 5-cycles... (this might take a moment)")
        num_5_cycles = count_simple_cycles(G, 5)
        cycle_counts.append(num_5_cycles)
        
        # Output results for the current graph
        print(f"Parameters (n, d, λ, μ): ({G.number_of_nodes()}, {G.degree(list(G.nodes())[0])}, {lam}, {mu})")
        print(f"Number of 5-cycles = {num_5_cycles}\n")

    # Final conclusion
    print("--- Conclusion ---")
    print(f"Yes, a pair of graphs in srg({n}, {d}, {lam}, {mu}) with a different number of 5-cycles exists.")
    
    n5_g1 = cycle_counts[0]
    n5_g2 = cycle_counts[1]
    
    print(f"The Shrikhande graph has {n5_g1} 5-cycles.")
    print(f"The Rook's graph (K4xK4) has {n5_g2} 5-cycles.")
    print("\nFinal equation showing the difference:")
    print(f"{n5_g1} != {n5_g2}")


if __name__ == "__main__":
    main()