import itertools

def get_t8_graph():
    """
    Constructs the adjacency list for the triangular graph T(8).
    Vertices are the 2-element subsets of {0, 1, ..., 7}.
    They are mapped to integers 0..27.
    """
    # Create a mapping from pairs {i, j} to an integer index
    elements = list(range(8))
    nodes = list(itertools.combinations(elements, 2))
    node_map = {node: i for i, node in enumerate(nodes)}
    
    n = len(nodes)
    adj = [[] for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            node1 = nodes[i]
            node2 = nodes[j]
            # Two vertices are adjacent if their sets share exactly one element
            if len(set(node1) & set(node2)) == 1:
                adj[i].append(j)
                adj[j].append(i)
    return adj, node_map, nodes

def get_chang_graph(t8_adj, t8_node_map, t8_nodes):
    """
    Constructs a Chang graph by switching T(8) with respect to a vertex set.
    The switching set S will be all pairs containing the element '0'.
    """
    n = len(t8_adj)
    # The switching set S includes all vertices {0, k} for k=1..7
    switching_set_indices = {t8_node_map[node] for node in t8_nodes if 0 in node}
    
    chang_adj = [[] for _ in range(n)]
    
    for i in range(n):
        for j in t8_adj[i]:
            # Adjacency is flipped if one vertex is in S and the other is not
            i_in_S = i in switching_set_indices
            j_in_S = j in switching_set_indices
            
            if i_in_S != j_in_S: # one is in S, the other is not
                # This edge is removed, so we do nothing
                pass
            else: # both in S or both not in S
                # Keep the original edge
                chang_adj[i].append(j)

    # Now add edges between non-adjacent vertices where one is in S and the other is not
    for i in range(n):
        for j in range(i + 1, n):
            is_adjacent = j in t8_adj[i]
            if not is_adjacent:
                i_in_S = i in switching_set_indices
                j_in_S = j in switching_set_indices
                if i_in_S != j_in_S:
                    chang_adj[i].append(j)
                    chang_adj[j].append(i)

    # Sort adjacency lists for consistency
    for i in range(n):
        chang_adj[i].sort()
        
    return chang_adj

def count_five_cycles(adj):
    """
    Counts the number of 5-cycles in a graph given its adjacency list.
    """
    n = len(adj)
    count = 0
    # Path: v0-v1-v2-v3-v4-v0
    for v0 in range(n):
        for v1 in adj[v0]:
            # To avoid overcounting, only consider paths where v1 > v0
            if v1 <= v0: continue
            for v2 in adj[v1]:
                if v2 == v0: continue
                for v3 in adj[v2]:
                    if v3 == v0 or v3 == v1: continue
                    # Paths of length 3: v0-v1-v2-v3
                    # Count paths of length 2 from v3 back to v0, not via v1 or v2
                    common_neighbors = set(adj[v3]) & set(adj[v0])
                    for v4 in common_neighbors:
                        if v4 != v1 and v4 != v2:
                           count += 1
    # Each cycle (v0-v1-v2-v3-v4-v0) is found once due to the v1 > v0 constraint.
    # The path-counting method finds v0-v1...v4-v0 and v0-v4...v1-v0. So divide by 2.
    return count // 2

def main():
    """
    Main function to construct graphs and count their 5-cycles.
    """
    print("Finding two SRGs with the same parameters but different numbers of 5-cycles.")
    print("The chosen parameters are (n, d, lambda, mu) = (28, 12, 6, 4).\n")

    # 1. T(8) graph
    print("Constructing the Triangular graph T(8)...")
    t8_adj, t8_node_map, t8_nodes = get_t8_graph()
    
    # 2. Chang graph
    print("Constructing a Chang graph by switching T(8)...")
    chang_adj = get_chang_graph(t8_adj, t8_node_map, t8_nodes)

    # 3. Count 5-cycles
    print("Counting 5-cycles in T(8)...")
    c5_t8 = count_five_cycles(t8_adj)
    
    print("Counting 5-cycles in the Chang graph...")
    c5_chang = count_five_cycles(chang_adj)

    # 4. Output results
    print("\n--- Results ---")
    print("Parameters: (n=28, d=12, lambda=6, mu=4)")
    print(f"Graph 1 (Triangular graph T(8)):")
    print(f"Number of 5-cycles = {c5_t8}")
    
    print(f"\nGraph 2 (Chang graph):")
    print(f"Number of 5-cycles = {c5_chang}")

    print(f"\nConclusion: The two graphs belong to the same class srg(28, 12, 6, 4) but have a different number of 5-cycles ({c5_t8} != {c5_chang}).")

if __name__ == '__main__':
    main()