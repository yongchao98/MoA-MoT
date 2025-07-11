import collections

def count_3_matchings(adj):
    """
    Counts the number of 3-matchings in a graph.

    The graph is represented by an adjacency list.
    This implementation runs in O(|V|*|E|) time, which is O(|V|^3) for dense graphs.
    Theoretically, faster subcubic algorithms exist (e.g., O(|V|^omega)),
    which is why k=3 is the answer to the problem.

    Args:
        adj (dict): An adjacency list representation of the graph.
                    Keys are vertices (0 to n-1), values are lists of neighbors.

    Returns:
        int: The number of 3-matchings in the graph.
    """
    n = len(adj)
    if n < 6:
        print("Graph has fewer than 6 vertices, so no 3-matching is possible.")
        print("0 / 3 = 0")
        return 0

    degrees = {i: len(adj.get(i, [])) for i in range(n)}
    m = sum(degrees.values()) // 2
    
    # Adjacency matrix for quick neighbor lookups
    adj_matrix = [[False] * n for _ in range(n)]
    edge_list = []
    for u in range(n):
        for v in adj.get(u, []):
            adj_matrix[u][v] = True
            if u < v:
                edge_list.append((u, v))

    total_2_matchings_in_subgraphs = 0

    # Iterate over each edge (u, v) once
    for u, v in edge_list:
        # Consider the subgraph G' = G - {u, v}
        
        # Number of edges in G'
        m_prime = m - degrees[u] - degrees[v] + 1
        
        # Calculate sum of C(deg_G'(w), 2) for all w in G'
        sum_binom_deg_prime = 0
        for w in range(n):
            if w == u or w == v:
                continue
            
            # Calculate degree of w in G'
            adj_to_uv = 0
            if adj_matrix[w][u]:
                adj_to_uv += 1
            if adj_matrix[w][v]:
                adj_to_uv += 1
            
            deg_prime_w = degrees[w] - adj_to_uv
            
            if deg_prime_w >= 2:
                sum_binom_deg_prime += deg_prime_w * (deg_prime_w - 1) // 2

        # Number of 2-matchings in G'
        # Formula: C(m', 2) - sum(C(deg'(w), 2) for w in G')
        num_2_matchings_in_g_prime = 0
        if m_prime >= 2:
            num_2_matchings_in_g_prime = (m_prime * (m_prime - 1) // 2) - sum_binom_deg_prime
        
        total_2_matchings_in_subgraphs += num_2_matchings_in_g_prime

    # Each 3-matching is counted 3 times (once for each of its edges)
    # The final equation is: 3 * (#3-matchings) = total_2_matchings_in_subgraphs
    num_3_matchings = total_2_matchings_in_subgraphs // 3
    
    print(f"The number of 3-matchings is derived from the equation:")
    print(f"{total_2_matchings_in_subgraphs} / 3 = {num_3_matchings}")
    
    return num_3_matchings

if __name__ == '__main__':
    # Example usage with a 6-cycle graph C6
    # Vertices 0, 1, 2, 3, 4, 5
    # Edges (0,1), (1,2), (2,3), (3,4), (4,5), (5,0)
    # This graph has two 3-matchings: {(0,1), (2,3), (4,5)} and {(1,2), (3,4), (5,0)}
    
    print("Counting 3-matchings in a 6-cycle graph (C6):")
    c6_adj = {
        0: [1, 5],
        1: [0, 2],
        2: [1, 3],
        3: [2, 4],
        4: [3, 5],
        5: [4, 0]
    }
    count_3_matchings(c6_adj)

    print("\n" + "="*30 + "\n")
    
    # Example with a complete graph K5 (no 3-matchings possible)
    print("Counting 3-matchings in a complete graph on 5 vertices (K5):")
    k5_adj = {i: [j for j in range(5) if i != j] for i in range(5)}
    count_3_matchings(k5_adj)