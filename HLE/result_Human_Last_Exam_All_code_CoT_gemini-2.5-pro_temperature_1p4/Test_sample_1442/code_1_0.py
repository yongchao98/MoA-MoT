import collections

def count_k_matchings(edges, k):
    """
    Counts the number of k-matchings in a graph represented by a list of edges.
    This is a recursive backtracking algorithm.
    """
    # Base case: if we need to find 0 more edges, we have found one valid matching.
    if k == 0:
        return 1
    
    # If we still need edges but have no more to choose from, this path is not a solution.
    if not edges:
        return 0

    # Recursive step:
    # 1. Take the first edge from the list.
    edge = edges[0]
    u, v = edge
    
    # 2. Count matchings that *do not* include this edge.
    #    We recurse on the rest of the edges.
    count_without_edge = count_k_matchings(edges[1:], k)
    
    # 3. Count matchings that *do* include this edge.
    #    We need to find k-1 more edges from the set of edges that are disjoint from the current one.
    remaining_edges_disjoint = [e for e in edges[1:] if u not in e and v not in e]
    count_with_edge = count_k_matchings(remaining_edges_disjoint, k - 1)
    
    return count_without_edge + count_with_edge

def get_edge_list(adj_list):
    """Converts an adjacency list to a unique list of edges."""
    edges = set()
    for u, neighbors in adj_list.items():
        for v in neighbors:
            # To ensure each edge is added once, we use a consistent order (u < v).
            if u < v:
                edges.add((u, v))
    return list(edges)

def solve():
    """
    Defines the two graphs, counts their 3-matchings, and prints the results.
    """
    # G1 is the graph of the cube, Q_3
    G1_adj = {
        0: [4, 5, 6], 1: [4, 5, 7], 2: [4, 6, 7], 3: [5, 6, 7],
        4: [0, 1, 2], 5: [0, 1, 3], 6: [0, 2, 3], 7: [1, 2, 3]
    }

    # G2 is the "twisted cube" graph
    G2_adj = {
        0: [5, 6, 7], 1: [4, 5, 7], 2: [4, 6, 7], 3: [4, 5, 6],
        4: [1, 2, 3], 5: [0, 1, 3], 6: [0, 2, 3], 7: [0, 1, 2]
    }
    
    # Both G1 and G2 are bipartite, 3-regular graphs on 8 vertices.

    g1_edges = get_edge_list(G1_adj)
    g2_edges = get_edge_list(G2_adj)

    # Number of 3-matchings for G1
    num_3_matchings_g1 = count_k_matchings(g1_edges, 3)
    
    # Number of 3-matchings for G2
    num_3_matchings_g2 = count_k_matchings(g2_edges, 3)

    print(f"Graph G1 is the Cube graph Q_3.")
    print(f"Graph G2 is the Twisted Cube graph.")
    print(f"Both are bipartite, 3-regular graphs on 8 vertices.")
    print(f"The number of 3-matchings in G1 is: {num_3_matchings_g1}")
    print(f"The number of 3-matchings in G2 is: {num_3_matchings_g2}")

    if num_3_matchings_g1 != num_3_matchings_g2:
        print("\nThe number of 3-matchings is not the same for both graphs.")
        print("Therefore, two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")
    else:
        print("\nThe number of 3-matchings is the same for both graphs in this case.")

solve()