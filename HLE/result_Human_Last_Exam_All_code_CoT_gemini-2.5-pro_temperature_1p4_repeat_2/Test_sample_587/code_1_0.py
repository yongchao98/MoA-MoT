import collections

def count_p4_subgraphs(graph):
    """
    Counts the number of subgraphs isomorphic to a path of 4 vertices (P_4) in a given graph.
    The graph is represented as an adjacency list.
    
    A P_4 subgraph is a simple path v1-v2-v3-v4 where all vertices are distinct.
    
    The number of injective homomorphisms (ordered paths) is counted first.
    An ordered path is specified by a sequence of 4 distinct vertices (v1, v2, v3, v4)
    such that (v1,v2), (v2,v3), and (v3,v4) are edges in the graph.
    
    The number of P_4 subgraphs is the number of these ordered paths divided by the
    size of the automorphism group of P_4, which is 2 (the identity and the reversal).
    """
    
    # Ensure graph is represented as adjacency sets for efficient lookups
    adj = {u: set(v) for u, v in graph.items()}
    
    # Initialize a default dictionary for nodes not in the primary keys
    adj = collections.defaultdict(set, adj)
    
    # Count injective homomorphisms (ordered v1-v2-v3-v4 paths)
    injective_hom_count = 0
    
    # The algorithm iterates over the middle edge (v2, v3) of the path.
    for v2 in adj:
        for v3 in adj[v2]:
            # Now we have the directed edge v2 -> v3
            
            # Count potential v1's. v1 must be a neighbor of v2, but not v3 itself.
            # (since path must be simple).
            count_v1 = len(adj[v2] - {v3})
            
            # Count potential v4's. v4 must be a neighbor of v3, but not v2 itself.
            count_v4 = len(adj[v3] - {v2})
            
            # These counts would be multiplied if we knew v1 and v4 are different.
            # We must subtract the cases where v1 == v4.
            # This happens if v1 (or v4) is a common neighbor of v2 and v3.
            common_neighbors = len(adj[v2].intersection(adj[v3]))
            
            # Number of ordered pairs (v1, v4) where v1 and v4 are distinct.
            term = count_v1 * count_v4 - common_neighbors
            
            injective_hom_count += term
            
    # The number of subgraphs is the homomorphism count divided by |Aut(P_4)| = 2.
    num_subgraphs = injective_hom_count // 2
    
    print(f"Input graph has {len(adj)} vertices and {sum(len(v) for v in adj.values()) // 2} edges.")
    print("The number of ordered 4-vertex paths (injective homomorphisms from P4) is {}.".format(injective_hom_count))
    print("The automorphism group of P4 has size 2.")
    print("Final equation:")
    print(f"Number of P4 subgraphs = {injective_hom_count} / 2 = {num_subgraphs}")

# Example usage with a complete graph on 5 vertices (K_5)
# In a K_n, any sequence of k distinct vertices forms a P_k.
# Number of ordered P4 paths in K_5 is P(5,4) = 5*4*3*2 = 120.
# Number of P4 subgraphs in K_5 = 120 / 2 = 60.
n = 5
k5_graph = {i: [j for j in range(n) if i != j] for i in range(n)}

count_p4_subgraphs(k5_graph)