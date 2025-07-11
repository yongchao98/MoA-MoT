import collections

def reduce_long_cycle_to_disjoint_cycles(g_prime, k_prime):
    """
    Demonstrates the W[1]-hardness of DisjointCycles by reducing Long Cycle to it.

    The Long Cycle problem (given G', k', find a simple cycle of length >= k')
    is a known W[1]-complete problem. We reduce it to DisjointCycles
    (given G, k, find k vertex-disjoint simple cycles of length >= k).

    Args:
        g_prime (dict): An adjacency list representation of the graph G'.
                        Keys are vertices (integers), values are lists of neighbors.
        k_prime (int): The parameter k' for the Long Cycle problem.
    """
    print("--- Reduction from Long Cycle to DisjointCycles ---")
    print(f"Input Long Cycle instance: (G', k')")
    print(f"k' = {k_prime}")
    # Using sorted keys for deterministic output
    print(f"G' adjacency list: {{_ for _ in { {k: sorted(v) for k, v in sorted(g_prime.items())} } }}")
    num_vertices_g_prime = len(g_prime)
    print(f"G' has {num_vertices_g_prime} vertices and {sum(len(v) for v in g_prime.values()) // 2} edges.\n")

    # The new parameter k for DisjointCycles is the same as k_prime.
    k = k_prime

    # The new graph G is the disjoint union of G' and (k-1) cycles of length k.
    g = collections.defaultdict(list)
    
    # Copy G' into G
    for u, neighbors in g_prime.items():
        g[u].extend(neighbors)

    # Add k-1 new cycles, each of length k
    # Start vertex labeling for the new cycles after the last vertex of G'
    if not g_prime:
        next_vertex = 0
    else:
        next_vertex = max(g_prime.keys()) + 1
        
    total_added_vertices = 0
    for i in range(k - 1):
        # Create a cycle of length k
        cycle_start_node = next_vertex
        if k > 0:
            for j in range(k):
                u = cycle_start_node + j
                v = cycle_start_node + ((j + 1) % k)
                g[u].append(v)
                g[v].append(u)
            next_vertex += k
            total_added_vertices += k

    print("Output DisjointCycles instance: (G, k)")
    print(f"k = {k}")
    # Using sorted keys for deterministic output
    print(f"G adjacency list: {{_ for _ in { {k: sorted(v) for k, v in sorted(g.items())} } }}")
    
    final_num_vertices = len(g)
    final_num_edges = sum(len(v) for v in g.values()) // 2
    
    print("\n--- Summary of Reduction ---")
    print(f"The constructed graph G is a disjoint union of the original graph G' and {k-1} new cycles of length {k}.")
    print("G has k disjoint cycles of length >= k if and only if G' has a cycle of length >= k'.")
    print(f"Number of vertices in G: {num_vertices_g_prime} (from G') + {k-1} * {k} (from new cycles) = {final_num_vertices}")
    print(f"This reduction shows that DisjointCycles is W[1]-hard, as Long Cycle is W[1]-complete.")
    

if __name__ == '__main__':
    # Example usage:
    # Let's take an instance of Long Cycle where G' is a diamond graph
    # and k' = 4. The graph has a 4-cycle (0-1-2-3-0).
    g_prime_example = {
        0: [1, 3],
        1: [0, 2, 3],
        2: [1, 3],
        3: [0, 1, 2]
    }
    k_prime_example = 4
    
    reduce_long_cycle_to_disjoint_cycles(g_prime_example, k_prime_example)

<<<B>>>