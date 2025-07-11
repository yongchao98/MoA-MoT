import networkx as nx

def long_cycle_to_disjoint_cycles_reduction(G_prime, k_prime):
    """
    This function demonstrates the reduction from an instance of LongCycle
    to an instance of DisjointCycles.

    Args:
        G_prime (nx.Graph): The input graph for the LongCycle problem.
        k_prime (int): The required cycle length for the LongCycle problem.

    Returns:
        A tuple (G, k) representing the instance for the DisjointCycles problem.
    """
    # The parameter k for DisjointCycles is the same as k_prime for LongCycle
    k = k_prime

    # Create the new graph G, starting with a copy of G_prime
    G = G_prime.copy()
    num_nodes_in_G_prime = G.number_of_nodes()

    # Add k-1 new, disjoint cycles of length k
    for i in range(k - 1):
        # Nodes for the new cycle. Offset by existing nodes and previous new cycles.
        start_node_idx = num_nodes_in_G_prime + i * k
        cycle_nodes = range(start_node_idx, start_node_idx + k)
        
        # Add edges to form a cycle
        for j in range(k):
            u = cycle_nodes[j]
            v = cycle_nodes[(j + 1) % k]
            G.add_edge(u, v)

    return G, k

def main():
    # Example usage of the reduction
    # Let's create an instance of LongCycle
    # G' is a 5-cycle, C5
    G_prime = nx.cycle_graph(5)
    k_prime = 5
    
    print(f"Original LongCycle instance:")
    print(f"  G' has {G_prime.number_of_nodes()} nodes and {G_prime.number_of_edges()} edges.")
    print(f"  k' = {k_prime}")
    print(f"  G' has a cycle of length >= {k_prime}: True (the graph itself is a C5)\n")
    
    # Perform the reduction
    G, k = long_cycle_to_disjoint_cycles_reduction(G_prime, k_prime)
    
    print(f"Resulting DisjointCycles instance:")
    print(f"  G has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    print(f"  k = {k}")

    # For this example, G is the disjoint union of C5 and (5-1)=4 copies of C5.
    # It contains 5 disjoint cycles of length 5.
    # The problem asks: Does G have at least k disjoint cycles of length at least k?
    # In our case: Does G have at least 5 disjoint cycles of length at least 5?
    # The answer is True. This demonstrates the YES-instance mapping.
    # The original C5 is one cycle.
    # We added 4 other C5s.
    # Total is 5 cycles of length 5.
    num_added_cycles = k - 1
    print(f"  The reduction added {num_added_cycles} disjoint cycles of length {k}.")
    print(f"  The resulting graph G is expected to have a solution because the original G' had one.")
    
if __name__ == "__main__":
    main()
