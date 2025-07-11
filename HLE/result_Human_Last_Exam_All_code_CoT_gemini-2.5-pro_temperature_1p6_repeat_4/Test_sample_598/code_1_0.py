import collections

def build_rooks_graph():
    """
    Builds the adjacency list for the 4x4 Rook's graph.
    Vertices are integers 0-15. Vertex v corresponds to cell (r,c)
    where r = v // 4 and c = v % 4.
    Two vertices are adjacent if they are in the same row or column.
    """
    adj = collections.defaultdict(list)
    n = 16
    for i in range(n):
        for j in range(i + 1, n):
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            if r1 == r2 or c1 == c2:
                adj[i].append(j)
                adj[j].append(i)
    # Convert to a regular dict of lists for printing and consistency
    return {k: sorted(v) for k, v in adj.items()}

def build_shrikhande_graph():
    """
    Builds the adjacency list for the Shrikhande graph.
    Vertices are integers 0-15 representing elements of Z_4 x Z_4.
    v -> (i,j) where i = v // 4, j = v % 4.
    Adjacency is defined by the difference being in a specific connection set.
    """
    adj = collections.defaultdict(list)
    n = 16
    # The connection set S for the Cayley graph construction
    S = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    for v_idx in range(n):
        r1, c1 = v_idx // 4, v_idx % 4
        for dr, dc in S:
            r2 = (r1 + dr) % 4
            c2 = (c1 + dc) % 4
            neighbor_idx = 4 * r2 + c2
            adj[v_idx].append(neighbor_idx)
    # Convert to a regular dict of sorted lists
    return {k: sorted(v) for k, v in adj.items()}

def count_cycles_of_length_5(adj):
    """
    Counts the number of simple 5-cycles in a graph.
    It iterates through all paths of length 4 (v0-v1-v2-v3-v4) and
    checks if an edge exists between v4 and v0 to close the cycle.
    It ensures all vertices in the cycle are distinct.
    The final count is divided by 10 to correct for overcounting
    (5 starting points and 2 directions for each cycle).
    """
    n = len(adj)
    count = 0
    for v0 in range(n):
        # Path: v0 -> v1
        for v1 in adj[v0]:
            # Path: v0 -> v1 -> v2
            for v2 in adj[v1]:
                if v2 == v0: continue
                # Path: v0 -> v1 -> v2 -> v3
                for v3 in adj[v2]:
                    if v3 == v1 or v3 == v0: continue
                    # Path: v0 -> v1 -> v2 -> v3 -> v4
                    for v4 in adj[v3]:
                        if v4 == v2 or v4 == v1: continue
                        # Check for a closing edge v4 -> v0
                        if v4 in adj[v0]:
                            # Ensure all vertices are distinct to form a simple cycle
                            if len({v0, v1, v2, v3, v4}) == 5:
                                count += 1
    
    # Each cycle is counted 10 times (5 starting vertices, 2 directions)
    return count // 10

def main():
    """
    Main function to construct the graphs, count their 5-cycles, and print the results.
    """
    # The parameters are (n, d, lambda, mu) = (16, 6, 2, 2)
    n = 16
    d = 6
    lambda_ = 2
    mu = 2

    # Create the two graphs
    rooks_graph = build_rooks_graph()
    shrikhande_graph = build_shrikhande_graph()
    
    # Count the 5-cycles in each graph
    rooks_5_cycles = count_cycles_of_length_5(rooks_graph)
    shrikhande_5_cycles = count_cycles_of_length_5(shrikhande_graph)
    
    print(f"Both graphs are in the class srg({n}, {d}, {lambda_}, {mu}).")
    print("-" * 40)
    print("Graph 1: 4x4 Rook's Graph")
    print(f"Number of 5-cycles: {rooks_5_cycles}")
    print("-" * 40)
    print("Graph 2: Shrikhande Graph")
    print(f"Number of 5-cycles: {shrikhande_5_cycles}")
    print("-" * 40)
    
    if rooks_5_cycles != shrikhande_5_cycles:
        print("Conclusion: The two graphs have a different number of 5-cycles.")
    else:
        print("Conclusion: The two graphs have the same number of 5-cycles.")


if __name__ == '__main__':
    main()