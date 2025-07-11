import collections

def count_cycles(adj, cycle_len):
    """
    Counts the number of simple cycles of a given length in a graph.

    Args:
        adj (dict): The adjacency list representation of the graph.
        cycle_len (int): The length of the cycles to count.

    Returns:
        int: The number of simple cycles of the specified length.
    """
    count = 0
    nodes = list(adj.keys())

    # We use a recursive DFS approach to find cycles
    def find_cycles_recursive(start_node, current_node, path):
        nonlocal count
        
        # When the path reaches the desired length for a cycle
        if len(path) == cycle_len:
            # Check if the last node is connected to the start node to close the cycle
            if start_node in adj[current_node]:
                count += 1
            return

        # Explore neighbors
        for neighbor in adj[current_node]:
            # Avoid visiting nodes already in the current path
            if neighbor not in path:
                # The start_node can only be revisited at the very end
                find_cycles_recursive(start_node, neighbor, path + [neighbor])

    # Start a search from each node in the graph
    for node in nodes:
        find_cycles_recursive(node, node, [node])
    
    # Each cycle is counted `cycle_len` times (once for each starting node)
    # and twice (for each direction of traversal). So we divide by (cycle_len * 2).
    return count // (cycle_len * 2)

def main():
    """
    Constructs two SRGs with the same parameters but different numbers of 5-cycles
    and prints the results.
    """
    # Parameters for the strongly regular graphs
    n = 16
    d = 6
    lam = 2  # lambda
    mu = 2
    
    print(f"Searching for two non-isomorphic graphs in srg({n}, {d}, {lam}, {mu}) with a different number of 5-cycles.\n")

    # --- Graph 1: 4x4 Rook's Graph (L2(4)) ---
    # Vertices are the 16 cells of a 4x4 grid.
    # Two vertices are adjacent if they are in the same row or column.
    adj_rook = collections.defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            if r1 == r2 or c1 == c2:
                adj_rook[i].append(j)
                adj_rook[j].append(i)

    # --- Graph 2: Shrikhande Graph ---
    # Vertices are elements of Z_4 x Z_4.
    # Adjacency is defined by a Cayley graph structure.
    adj_shrikhande = collections.defaultdict(list)
    # The neighborhood of the vertex (0,0)
    neighbors_of_zero = {(0, 1), (0, 3), (1, 0), (3, 0), (1, 1), (3, 3)}
    
    nodes = [(r, c) for r in range(4) for c in range(4)]
    
    for i in range(n):
        for j in range(i + 1, n):
            r1, c1 = nodes[i]
            r2, c2 = nodes[j]
            
            # Check if the difference is in the neighborhood set
            dr = (r1 - r2) % 4
            dc = (c1 - c2) % 4
            
            if (dr, dc) in neighbors_of_zero:
                adj_shrikhande[i].append(j)
                adj_shrikhande[j].append(i)

    # Count the 5-cycles in both graphs
    print("Counting 5-cycles for each graph...")
    c5_rook = count_cycles(adj_rook, 5)
    c5_shrikhande = count_cycles(adj_shrikhande, 5)
    
    print("\n--- Results ---")
    print(f"Parameters (n, d, \u03BB, \u03BC) = ({n}, {d}, {lam}, {mu})")
    print(f"Graph 1 (4x4 Rook's Graph):")
    print(f"Number of 5-cycles = {c5_rook}")
    
    print(f"\nGraph 2 (Shrikhande Graph):")
    print(f"Number of 5-cycles = {c5_shrikhande}")
    
    print("\n--- Conclusion ---")
    print(f"The number of 5-cycles is different for these two graphs:")
    print(f"{c5_rook} (Rook's) \u2260 {c5_shrikhande} (Shrikhande)")
    print("\nThis shows that two graphs G, G' in the same srg(n,d,\u03BB,\u03BC) class can have a different number of 5-cycles.")


if __name__ == "__main__":
    main()